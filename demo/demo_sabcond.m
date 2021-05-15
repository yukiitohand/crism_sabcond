% demo script for performing atmospheric correction and de-noising
crism_init;

% Enter the list of observation IDs.
obs_id_list = {'9036'};

%%
% Download images and calibration.
for i=1:length(obs_id_list)
    obs_id = obs_id_list{i};
    crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L','Download_DDR',2,...
       'Download_EDRSCDF',2,'Download_TRRIF',2,'Download_TRRRA',2);
       % 'Download_TER',2);
end

for i=1:length(obs_id_list)
    obs_id = obs_id_list{i};
    crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');
    TRRIFdata = CRISMdata(crism_obs.info.basenameIF,crism_obs.info.dir_trdr);
    TRRIFdata.load_basenamesCDR('Dwld',2);
    crism_calibration_IR_v2(obs_id,'save_memory',true,'mode','yuki2', ...
        'version','B','skip_ifexist',0,'force',1);
end

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
save_pdir = './resu/';
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
% storage_saving_level = 'NORMAL';
% additional_suffix  = '';
% interleave_out     = 'lsb';
% interleave_default = 'lsb';
% subset_columns_out = false;
% Alib_out           = false;
% do_crop_bands      = false;

%  ## INPUT IMAGE OPTIONS #------------------------------------------------
opt_img      = 'TRRB';
% img_cube     = [];
% img_cube_band_inverse = [];
% dir_yuk      = crism_env_vars.dir_YUK; % TRRY_PDIR
% ffc_counter  = 1;

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
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


%% Main routine

for i=1:length(obs_id_list)
    obs_id = obs_id_list{i};
    result = ...
        sabcondv3_pub_water_ice_test(obs_id,3,...'t_mode',t_mode,'lambda_a',lambda_a,...
            'opt_img',opt_img,'OPTBP',optBP,'nIter',nIter,'Bands_Opt',bands_opt,...
            'PROC_MODE',proc_mode,'precision',precision);
    fprintf('water_ice_result: %d\n',result.presence_H2Oice);
    fprintf('water_ice_exist: %f\n',mean(result.presence_H2Oice_columns,'omitnan'));
    if result.presence_H2Oice
        sabcondv3_pub(obs_id,'t_mode',t_mode,'lambda_a',lambda_a,'opt_img',opt_img,...
            'OPTBP',optBP,'nIter',nIter,'Bands_Opt',bands_opt,'SAVE_PDIR',save_pdir,...
            'additional_suffix',additional_suffix,'force',force,'skip_ifexist',skip_ifexist,...
            'PROC_MODE',proc_mode,'precision',precision,'OPT_ICELIB',3);
    else
       sabcondv3_pub(obs_id,'t_mode',t_mode,'lambda_a',lambda_a,'opt_img',opt_img,...
           'OPTBP',optBP,'nIter',nIter,'Bands_Opt',bands_opt,'SAVE_PDIR',save_pdir,...
           'additional_suffix',additional_suffix,'force',force,'skip_ifexist',skip_ifexist,...
           'PROC_MODE',proc_mode,'precision',precision,...
           'OPT_CRISMTYPELIB',4,'OPT_SPLIBUSGS',6);
    end
    
%     sabcondv3_pub(obs_id,'t_mode',2,'lambda_a',0.01,'opt_img','TRRB',...
%        'OPTBP','pri','nIter',5,'Bands_Opt',bands_opt,'SAVE_PDIR',save_pdir,...
%        'additional_suffix','woice','force',force,'skip_ifexist',skip_ifexist,...
%        'PROC_MODE',proc_mode,'precision',precision,...
%        'OPT_CRISMTYPELIB',4,'OPT_SPLIBUSGS',6);
    
    
end

%%
