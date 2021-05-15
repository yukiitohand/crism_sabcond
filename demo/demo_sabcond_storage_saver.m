% demo script for performing atmospheric correction and de-noising
crism_init;

% Enter the list of observation IDs.
obs_id_list = {'9036'};

%%
% OPTIONS for sabcond
%
% -------------------------------------------------------------------------
%                               Major options
% -------------------------------------------------------------------------
%
% ## PROCESSING OPTIONS #--------------------------------------------------
precision  = 'single'; 
% {'single','double'}
% in my test, single was much faster and result is quite similar.

proc_mode  = 'GPU_BATCH_2'; 
% {'CPU_1','GPU_1','CPU_2','GPU_2','GPU_BATCH_2'}
% 'GPU_BATCH_2' is the fastest mode, if GPU is not used by others.
% Slightly different algorithms between {'CPU_1','GPU_1'} and 
% {'CPU_2','GPU_2','GPU_BATCH_2'}. Latter is supposed to be better.

batch_size = 10; 
% column size for which processing is performed. Valid only if 
% 'GPU_BATCH_*' mode is selected.

% ## I/O OPTIONS #---------------------------------------------------------
storage_saving_level = 'Highest';
% string {'Highest','Normal'}
% determine how much to save storage, no effect under save_file=1
%    Normal  - All the byproducts are saved
%    Highest - Only nr_ds and mdl_ds are saved, bands are also scropped.
% (default) Normal

save_pdir = '';
% character, string
% root directory path where the processed data are stored. The processed 
% image will be saved at <SAVE_PDIR>/CCCNNNNNNNN, where CCC the class type 
% of the obervation and NNNNNNNN is the observation id. It doesn't matter 
% if trailing slash is there or not.

save_dir_yyyy_doy = 0;
% Boolean
% if true, processed images are saved at <SAVE_PDIR>/YYYY_DOY/CCCNNNNNNNN,
% otherwise, <SAVE_PDIR>/CCCNNNNNNNN.

force = 0;
% Boolean
% if true, processing is forcefully performed and all the existing images 
% will overwritten. Otherwise, you will see a prompt asking whether or not 
% to continue and overwrite images or not when there alreadly exist 
% processed images.

skip_ifexist = 1;
% Boolean
% if true, processing will be automatically skipped if there already exist 
% processed images. No prompt asking whether or not to continue and 
% overwrite images or not.

additional_suffix = '';
% character, string
% any additional suffix added to the name of processd images.

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
bands_opt    = 4;
% integer
% Magic number for wavelength channels to be used. This is the input for 
% [bands] = genBands(bands_opt)

%%
% -------------------------------------------------------------------------
%                               Minor Options
% -------------------------------------------------------------------------

% ## I/O OPTIONS #---------------------------------------------------------
% save_file          = true;
% additional_suffix  = '';
% interleave_out     = 'lsb';
% interleave_default = 'lsb';
% subset_columns_out = false;
% Alib_out           = false;

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
opt_img      = 'TRRB';
% dir_yuk      = crism_env_vars.dir_YUK; % TRRX_PDIR
% ffc_counter  = 1;
% line_idxes   = [];                     % LINES
% column_idxes = [];
% mt           = 'sabcondpub_v1';        % METHODTYPE
optBP        = 'pri';                  %{'pri','all','none'}
% verbose      = 0;
% column_skip  = 1;
% weight_mode  = 0;
% lambda_update_rule = 'L1SUM';
% th_badspc    = 0.8;
% 
% ## TRANSMISSION SPECTRUM OPTIONS #---------------------------------------
t_mode = 2;
% obs_id_T = '';
% varargin_T = {};

% ## LIBRARY OPTIONS #-----------------------------------------------------
% cntRmvl         = 1;
% optInterpid     = 1;

% following are partially controlled by automatic detection of water ice,
% change them with caution if you want.
% optCRISMspclib  = 1;
% optRELAB        = 1;
% optUSGSsplib    = 1;
% optCRISMTypeLib = 2;
% opticelib       = '';
% 
% ## SABCONDC OPTIONS #----------------------------------------------------
nIter = 5;
lambda_a = 0.01;


%%
% Download images and calibration.
% You need to include all the processing before going to a next image.
for i=1:length(obs_id_list)
    obs_id = obs_id_list{i};
    
    %---------------------------------------------------------------------%
    % Download images
    %---------------------------------------------------------------------%
    % First download images and labels from the observation folder.
    % EDR SC DF images, DDR images, LBL of TRR3 I/F RA, and TRR3 HKP files
    % are downloaded. 
    % Latest downloader use a cached index.html in the local storage to get
    % the list of files. If the remote server was configured different, it
    % might not be able to find the list of files. Setting
    %   DWLD_INDEX_CACHE_UPDATE=1
    % always update index.html, so it is safe to apply this option. This is
    % more critical for downloading CDR...
    crism_obs_info = crism_get_obs_info(obs_id,'SENSOR_ID','L', ...
        'Download_DDR',2,'Download_EDRSCDF',2, ...
        'Download_TRRIF',2,'ext_trrif','lbl', ...
        'Download_TRRRA',2,'ext_trrra','lbl','Download_TRRRAHKP',2,...
        'DWLD_INDEX_CACHE_UPDATE',true);

    % Download CDR data necessary for calibration and later processing.
    % fnamesCDR_local gets the filenames of all the relevant CDR files
    % existing in the local storage (used for deleting once the processing
    % is done)
    TRRIFdata = CRISMdata(crism_obs.info.basenameIF,crism_obs.info.dir_trdr);
    [fnamesCDR_local] = TRRIFdata.load_basenamesCDR('Dwld',2,'INDEX_CACHE_UPDATE',true);
    
    %---------------------------------------------------------------------%
    % Calibration
    %---------------------------------------------------------------------%
    % Download=0, no additional downloading is performed. Currently this is
    % supported with mode=yuki2. yuki3 and yuki4 may need additional data
    % to be downloaded.
    [Yif] = crism_calibration_IR_v2(obs_id,'save_memory',true,'mode','yuki2', ...
        'version','B','skip_ifexist',0,'force',1, 'Dwld',0);
    
    %---------------------------------------------------------------------%
    % Atmospheric correction and denoising.
    %---------------------------------------------------------------------%
    result = ...
        sabcondv3_pub_water_ice_test(obs_id,3,...'t_mode',t_mode,'lambda_a',lambda_a,...
            'opt_img',opt_img,'img_cube',Yif,'img_cube_band_inverse',0, ...
            'OPTBP',optBP,'nIter',nIter,'Bands_Opt',bands_opt,...
            'PROC_MODE',proc_mode,'precision',precision);
    fprintf('water_ice_result: %d\n',result.presence_H2Oice);
    fprintf('water_ice_exist: %f\n',mean(result.presence_H2Oice_columns,'omitnan'));
    if result.presence_H2Oice
        sabcondv3_pub(obs_id,'t_mode',t_mode,'lambda_a',lambda_a, ...
            'opt_img',opt_img,'img_cube',Yif,'img_cube_band_inverse',0, ...
            'OPTBP',optBP,'nIter',nIter,'Bands_Opt',bands_opt,'SAVE_PDIR',save_pdir,...
            'additional_suffix',additional_suffix,'force',force,'skip_ifexist',skip_ifexist,...
            'PROC_MODE',proc_mode,'precision',precision,'OPT_ICELIB',3,...
            'STORAGE_SAVING_LEVEL',storage_saving_level);
    else
       sabcondv3_pub(obs_id,'t_mode',t_mode,'lambda_a',lambda_a, ...
           'opt_img',opt_img,'img_cube',Yif,'img_cube_band_inverse',0, ...
           'OPTBP',optBP,'nIter',nIter,'Bands_Opt',bands_opt,'SAVE_PDIR',save_pdir,...
           'additional_suffix',additional_suffix,'force',force,'skip_ifexist',skip_ifexist,...
           'PROC_MODE',proc_mode,'precision',precision,...
           'OPT_CRISMTYPELIB',4,'OPT_SPLIBUSGS',6, ...
           'STORAGE_SAVING_LEVEL',storage_saving_level);
    end
    
    %---------------------------------------------------------------------%
    % Cleaning. Delete files from the local storage.
    %---------------------------------------------------------------------%
    func_delete = @(dirpath,fnamecell) cellfun(@(x) delete(joinPath(dirpath,x)),fnamecell);
    func_delete(crism_obs_info.dir_trdr,crism_obs_info.fnameTRRwext_local);
    func_delete(crism_obs_info.dir_edr ,crism_obs_info.fnameEDRwext_local);
    func_delete(crism_obs_info.dir_ddr ,crism_obs_info.fnameDDRwext_local);
    
    % CDR
    % only CDRs updated daily (stored yyyy_doy folder) are removed.
    cdr_acros = fieldnames(fnamesCDR_local);
    for ii=1:length(cdr_acros)
        acro = cdr_acros{ii};
        fldtype = crism_assessCDRForderType(acro);
        switch fldtype
            case 1
                if iscell(TRRIFdata.dir_cdr.(acro))
                    for iii=1:length(fnamesCDR_local.(acro))
                        dirpath_cdr = TRRIFdata.dir_cdr.(acro){iii};
                        fnames_cdr  = fnamesCDR_local.(acro){iii};
                        func_delete(dirpath_cdr,fnames_cdr);
                    end
                elseif ischar(TRRIFdata.dir_cdr.(acro))
                    dirpath_cdr = TRRIFdata.dir_cdr.(acro);
                    fnames_cdr  = fnamesCDR_local.(acro);
                    func_delete(dirpath_cdr,fnames_cdr);
                end
        end
    end
    %---------------------------------------------------------------------%
end

%%
