function [out] = sabcondv3_pub(obs_id,varargin)
% [out] = sabcondv3_pub(obs_id,varargin)
%   This function performs simultaneous de-noising and atmospheric
%   correction of CRISM data. For any technical information, refer the
%   following paper:
%
%       A new atmospheric correction and de-noising method for CRISM data
%       Y. Itoh and M. Parente
%       Submitted. Pre-print available. 20 pages.
%       http://rhogroup.org/images/downloads/crism-atmcorr-denoising.pdf
%
%   By default, three kinds of images are saved in the ENVI format:
%
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX].img
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX]_nr.img
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX]_ori.img
%
%   These are accompanied with header files. The first image is the base
%   product where bad entires are substituted with NaNs, the "*_nr" image 
%   is the corrected image whose bad entries are substituted with our best 
%   guess, and the "*_ori" image is the corrected image whose bad entries 
%   are left. [BASENAME_INPUT] depends on the option 'OPT_IMG'. You can 
%   also customize the output image names by specifying the options 
%   'METHODTYPE' and 'ADDITIONAL_SUFFIX'. 
%
%   All the OPTIONAL_parameter settings are saved in
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX]_settings.txt
%   and ancillary information is saved 
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX]_AB.img
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX]_Bg.img
%       [BASENAME_INPUT]_[METHODTYPE]_[ADDITIONAL_SUFFIX].mat
%   "*_AB" image is the model absorption image, "*_Bg" is the estimated 
%   background image, ".mat" stores other ancillary information.
%
%   Processed images will be saved at a specified folder. Read I/O options
%   in Optional Parameters later for how to specify folders. You can also
%   select an option not to save the processed images (See 'SAVE_FILE')
%
% INPUT Parameters
%   obs_id : string, 
%       observation ID of the image to be processed. Currently, EPF 
%       measurements accompanied with scene measurements are not supported.
%
% OUTPUT Parameters
%   out: struct, 
%       storing processed images and byproducts. 
%           out.Yif_cor     = Yif_cor;
%           out.Yif_cor_nr  = Yif_cor_nr;
%           out.Yif_cor_ori = Yif_cor_ori;
%           out.AB_est      = AB_est;
%           out.Bg_est      = Bg_est;
%           out.T_est       = T_est;
%           out.Yif_isnan    = Yif_isnan;
%           out.Valid_pixels = Valid_pixels;
%           out.WA          = WA;
%           out.lines       = line_idxes;
%           out.columns     = Columns;
%           out.bands       = bands;
%           out.interleave_out     = interleave_out;
%           out.subset_columns_out = subset_columns_out;
%       **ADDITIONAL FIELD**
%        For PROC_MODE either of {'CPU_1','GPU_1'}
%           out.ancillaries  = ancillaries;
%        For PROC_MODE either of {'GPU_BATCH_1'}
%           
%           out.X            = X;
%
% OPTIONAL Parameters 
%  ## I/O OPTIONS  #-------------------------------------------------------
%   'SAVE_FILE': boolean
%       whether or not to save processed images. If true, two optioal 
%       parameters 'FORCE','SKIP_IFEXIST' have no effect.
%       (default) true
%   'SAVE_PDIR': any string
%       root directory path where the processed data are stored. The
%       processed image will be saved at <SAVE_PDIR>/CCCNNNNNNNN, where CCC
%       the class type of the obervation and NNNNNNNN is the observation id.
%       It doesn't matter if trailing slash is there or not.
%       (default) './res/'
%   'SAVE_DIR_YYYY_DOY': boolean
%       if true, processed images are saved at 
%           <SAVE_PDIR>/YYYY_DOY/CCCNNNNNNNN,
%       otherwise, 
%           <SAVE_PDIR>/CCCNNNNNNNN.
%       (default) false
%   'FORCE': boolean
%       if true, processing is forcefully performed and all the existing
%       images will overwritten. Otherwise, you will see a prompt asking
%       whether or not to continue and overwrite images or not when there
%       alreadly exist processed images.
%       (default) false
%   'SKIP_IF_EXIST': boolean
%       if true, processing will be automatically skipped if there already 
%       exist processed images. No prompt asking whether or not to continue
%       and overwrite images or not.
%       (default) false
%   'ADDITIONAL_SUFFIX': any string,
%       any additional suffix added to the name of processd images.
%       (default) 'v1'
%   'INTERLEAVE_OUT': string, {'lsb','bls'}
%       interleave option of the images in the output parameter, out. This
%       is not the interleave used for saving processed images. 
%       'lsb': Line-Sample-Band
%       'bls': Band-Line-Sample
%       (default) 'lsb'
%   'SUBSET_COLUMNS_OUT': boolean,
%       whether or not to take only the images with the subset specified by
%       the optional parameter 'COLUMNS'. This is not the interleave used 
%       for saving processed images. 
%       (default) false
%   'ALIB_OUT': boolean,
%       whether or not to out to include the libraries used for the
%       processing.
%       (default) false
%
%  ## GENERAL SABCOND OPTIONS #--------------------------------------------
%   'OPT_IMG': string, {'IF','RA_IF','TRRY','TRRB','TRRC'}
%       type of input image to be used
%       (default) 'TRRB'
%   'TRRY_PDIR': any string
%       root directory path where {'TRRY','TRRB','TRRC'} images are stored.
%       If you choose either {'IF', 'RA_IF'} for 'OPT_IMG', you do not need
%       to specify this.
%   'FFC_IF_COUNTER': integer, {1,3}
%       Observation counter of the image to be processed. Used only for 
%       processing FFC images. 
%       (default) 1
%   'BANDS_OPT' : integer, {4}
%       an option for wavelength channels to use. This is the 
%       input for genBands()
%       (default) 4
%   'LINES': array
%       array, defining lines to be used for processing. If empty, then all
%       the lines are used.
%       (default) []
%   {'SAMPLES','COLUMNS','CLIST'}
%       array, defining columns to be processed. If empty, then all the
%       columns are processed.
%       (default) []
%   'METHODTYPE' : any string
%       type of the method
%       (default) sabcondv3_pub
%   'OPTBP' : string {'pri','all','none'}
%       option defining bad pixel information to be used.
%       (default) 'pri'
%   'VERBOSE': integer, {0,1,2}
%       option for how much information is displayed.
%
%  ## TRANSMISSION SPECTRUM OPTIONS #--------------------------------------
%    'T_MODE': integer
%       option for what kind of transmission spectrum frame is used.
%       (default) 2
%
%  ## LIBRARY OPTIONS #----------------------------------------------------
%   'CNTRMVL': boolean
%       whether or not to perform continuum removal on library spectra
%       (default) 1
%   'OPTINTERPID': {1}
%       option defining how the interpolation of invalid channels of the 
%       library spectra are dealt with.
%       (default) 1
%   'OPT_CRISMSPCLIB': integer, {1}
%       library option of the CRISM spectral library.
%       (default) 1
%   'OPT_RELAB': integer
%       integer, library option of RELAB
%       (default) 1
%   'OPT_SPLIBUSGS': integer
%       integer, library option of USGS splib
%       (default) 1
%   'OPT_CRISMTYPELIB': integer
%       integer, library option of CRISM MICA library
%       (default) 1
%   'OPT_ICELIB': integer
%       integer, library option of ICE library
%       (default) '' (no ice)
%
%  ## SABCONDC OPTIONS #---------------------------------------------------
%   'NITER': integer
%       the number of outer iteration
%       (default) 5
%   'LAMBDA_A': double array or scalar
%       trade-off parameter of the cost function.
%       (default) 0.01
%
%  ## PROCESSING OPTIONS #-------------------------------------------------
%   'PRECISION': string, {'single','double'}
%       percision with which processing is performed.
%       (default) 'double'
%   'PROC_MODE': string, {'CPU_1','GPU_1','CPU_2','GPU_2','GPU_BATCH_2'}
%       option for which processor is used CPU or GPU, and which algorithm
%       is used, formatted as PPP[_BATCH]_A, where PPP represents 
%       processor, and A is the index for algorithm. BATCH indicates 
%       (default) 'CPU_1'
%       **ALGORITHM description**
%           1: Bad pixel detected during processing are substitueded by
%              model values.
%           2: Bad pixel detected during processing are exactly ignored
%              during the processing.
%       Currently, multiple GPU processing is not supported.
%   'BATCH_SIZE': integer
%       column size for which processing is performed. Valid only if
%       'GPU_BATCH_*' mode is selected.
%       (default) 10

global crism_env_vars

% ## I/O OPTIONS #---------------------------------------------------------
save_file          = true;
save_pdir          = './resu/';
save_dir_yyyy_doy  = false;
force              = false;
skip_ifexist       = false;
additional_suffix  = 'v1';
interleave_out     = 'lsb';
interleave_default = 'lsb';
subset_columns_out = false;
Alib_out           = false;

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
opt_img      = 'TRRB';
dir_yuk      = crism_env_vars.dir_YUK; % TRRY_PDIR
ffc_counter  = 1;
bands_opt    = 4;
line_idxes   = [];                     % LINES
column_idxes = [];
mt           = 'sabcondpub_v1';        % METHODTYPE
optBP        = 'pri';                  %{'pri','all','none'}
verbose      = 0;

% ## TRANSMISSION SPECTRUM OPTIONS #---------------------------------------
t_mode = 2;

% ## LIBRARY OPTIONS #-----------------------------------------------------
cntRmvl         = 1;
optInterpid     = 1;
optCRISMspclib  = 1;
optRELAB        = 1;
optUSGSsplib    = 1;
optCRISMTypeLib = 2;
opticelib       = '';

% ## SABCONDC OPTIONS #----------------------------------------------------
nIter = 5;
lambda_a = 0.01;

% ## PROCESSING OPTIONS #--------------------------------------------------
precision  = 'double';
PROC_MODE  = 'CPU_1';
batch_size = 10;

% ## ETCETERA #------------------------------------------------------------
% gausssigma = 0.6;
% fltsize  = 5;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## I/O OPTIONS #---------------------------------------------
            case 'SAVE_FILE'
                save_file = varargin{i+1};
            case 'SAVE_PDIR'
                save_pdir = varargin{i+1};
            case 'SAVE_DIR_YYYY_DOY'
                save_dir_yyyy_doy = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case 'SKIP_IFEXIST'
                skip_ifexist = varargin{i+1};
            case 'ADDITIONAL_SUFFIX'
                additional_suffix = varargin{i+1};
            case 'INTERLEAVE_OUT'
                interleave_out = varargin{i+1};
            case 'SUBSET_COLUMNS_OUT'
                subset_columns_out = varargin{i+1};
            case 'ALIB_OUT'
                Alib_out = varargin{i+1};
                
            % ## GENERAL SABCOND OPTIONS #---------------------------------
            case 'OPT_IMG'
                opt_img = varargin{i+1};
            case 'TRRY_PDIR'
                dir_yuk = varargin{i+1};
            case 'FFC_IF_COUNTER'
                ffc_counter = varargin{i+1};
            case 'BANDS_OPT'
                bands_opt = varargin{i+1};
            case {'LINES','LLS'}
                line_idxes = varargin{i+1};
            case {'SAMPLES','COLUMNS','CLIST'}
                column_idxes = varargin{i+1};
            case 'METHODTYPE'
                mt = varargin{i+1};
            case 'OPTBP'
                optBP = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
                
            % ## TRANSMISSION SPECTRUM OPTIONS #---------------------------
            case 'T_MODE'
                t_mode = varargin{i+1};
                
            % ## LIBRARY OPTIONS #-----------------------------------------
            case 'CNTRMVL'
                cntRmvl = varargin{i+1};
            case 'OPTINTERPID'
                optInterpid = varargin{i+1};
            case 'OPT_CRISMSPCLIB'
                optCRISMspclib = varargin{i+1};
            case 'OPT_RELAB'
                optRELAB= varargin{i+1};
            case 'OPT_SPLIBUSGS'
                optUSGSsplib = varargin{i+1};
            case 'OPT_CRISMTYPELIB'
                optCRISMTypeLib = varargin{i+1};
            case 'OPT_ICELIB'
                opticelib = varargin{i+1};
                
            % ## SABCONDC OPTIONS #----------------------------------------
            case 'NITER'
                nIter = varargin{i+1};
            case 'LAMBDA_A'
                lambda_a = varargin{i+1};
                
            % ## PROCESSING OPTIONS #--------------------------------------
            case 'PRECISION'
                precision = varargin{i+1};
            case 'PROC_MODE'
                PROC_MODE = varargin{i+1};
            case 'BATCH_SIZE'
                batch_size = varargin{i+1};
                
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if save_file && force && skip_ifexist
    error('You are forcing or skipping? Not sure what you want');
end

if save_file && ~exist(save_pdir,'dir'), mkdir(save_pdir); end

bands = genBands(bands_opt);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
% libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,'');

switch upper(PROC_MODE)
    case {'CPU_1','CPU_2'}
        gpu = false;
    case {'GPU_1','GPU_2','GPU_BATCH_2'}
        gpu = true;
    otherwise
        error('Undefined PROC_MODE=%s',PROC_MODE);
end

% Check GPU ---------------------------------------------------------------
if gpu
    if gpuDeviceCount==0
        error('No GPU is detected. Use CPU option');
    elseif gpuDeviceCount > 1
        fprintf('Multiple GPUs are detected. The first one is used for processing');  
    end
end

%% Read image and ancillary data and format them for processing
crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');
switch upper(crism_obs.info.obs_classType)
    case {'FRT','HRL','HRS','FRS','ATO','MSP','HSP'}
        TRRIFdata = get_CRISMdata(crism_obs.info.basenameIF,'');
        TRRRAdata = get_CRISMdata(crism_obs.info.basenameRA,'');
    case {'FFC'}
        [TRRIFdata] = get_scene_CRISMdata_FFC(crism_obs.info.basenameIF,'',ffc_counter);
        [TRRRAdata] = get_scene_CRISMdata_FFC(crism_obs.info.basenameRA,'',ffc_counter);
end
TRRIF_is_empty = isempty(TRRIFdata);
if TRRIF_is_empty
    TRRIFdata = TRRRAdata;
end

[DFdata1,DFdata2] = get_DFdata4SC(TRRIFdata,crism_obs);

%%
%-------------------------------------------------------------------------%
% cheking the file exist or not.
%-------------------------------------------------------------------------%
if save_dir_yyyy_doy
    save_dir = joinPath(save_pdir,crism_obs.info.yyyy_doy,crism_obs.info.dirname);
else
    save_dir = joinPath(save_pdir,crism_obs.info.dirname);
end

if isempty(additional_suffix)
    suffix = mt;
else
    suffix = [mt '_' additional_suffix];
end

fprintf('suffix will be \n"%s"\n',suffix);

switch upper(opt_img)
    case 'IF'
        basename_cr = [crism_obs.info.basenameIF suffix];
    case 'RA_IF'
        basename_cr = [crism_obs.info.basenameRA '_IF' suffix];
    case {'TRRY','TRRB','TRRC'}
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = opt_img(4);
        basenameTRRY = get_basenameOBS_fromProp(prop);
        basename_cr = [basenameTRRY suffix];
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

fpath_cr = joinPath(save_dir,[basename_cr,'.img']);
if save_file && exist(fpath_cr,'file')
    if skip_ifexist
        return;
    elseif ~force
        flg = 1;
        while flg
            prompt = sprintf('There exists the image %s\n Do you want to continue to process and overwrite?(y/n)',fpath_cr);
            ow = input(prompt,'s');
            if any(strcmpi(ow,{'y','n'}))
                flg=0;
            else
                fprintf('Input %s is not valid.\n',ow);
            end
        end
        if strcmpi(ow,'n')
            fprintf('Process aborted...\n');
            diary off;
            return;
        elseif strcmpi(ow,'y')
            fprintf('processing continues and will overwrite...\n');
        end
    end
end


if save_file && ~exist(save_dir,'dir'), mkdir(save_dir); end

% open log file
username = char(java.lang.System.getProperty('user.name'));
fname = sprintf('log_%s_%s.txt',username,datetime('now','TimeZone','local','Format','yyyyMMdd'));
if save_file, diary(joinPath(save_dir,fname)); end

%% Read image and ancillary data and format them for processing
TRRIFdata.load_basenamesCDR();
WAdata = TRRIFdata.readCDR('WA'); WAdata.readimgi();
nLall = TRRIFdata.hdr.lines; nCall = TRRIFdata.hdr.samples; nBall = TRRIFdata.hdr.bands;
if isempty(line_idxes), line_idxes = 1:nLall; end
nL = length(line_idxes); nB = length(bands);
lBool = false(nLall,1); lBool(line_idxes) = true;
bBool = false(nBall,1); bBool(bands) = true;

basenameWA = WAdata.basename;
WAb = squeeze(WAdata.img(:,:,bands))';

switch upper(opt_img)
    case 'IF'
        if TRRIF_is_empty
            error('TRR I/F does not exist.');
        else
            Yif = TRRIFdata.readimgi();
        end
    case 'RA_IF'
        TRRRAIFdata = CRISMdataCAT([crism_obs.info.basenameRA '_IF'],crism_obs.info.dir_trdr,'ra_if');
        Yif = TRRRAIFdata.readimgi();
    case {'TRRY','TRRB','TRRC'}
        d_IoF = joinPath(dir_yuk, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        TRRYIFdata = CRISMdata(basenameTRRY,d_IoF);
        Yif = TRRYIFdata.readimgi();
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

Yif = Yif(line_idxes,:,bands);
Yif(Yif<=1e-8) = nan;
logYif = log(Yif);
logYif = permute(logYif,[3,1,2]);

fprintf('finish loading Image\n');

%%
% read bad pixel
[BPdata1,BPdata2,BPdata_post] = load_BPdataSC_fromDF(...
    TRRIFdata,DFdata1.basename,DFdata2.basename);
switch lower(optBP)
    case 'pri'
        [BP_pri1nan] = formatBPpri1nan(BPdata1,BPdata2,'band_inverse',true);
        [GP_pri] = convertBP1nan2GP1nan(BP_pri1nan);
        GP = permute(GP_pri(bands,:,:),[1,3,2]);
    case 'all'
        [BP_post1nan] = formatBP1nan(BPdata_post,'band_inverse',true);
        [GP_all] = convertBP1nan2GP1nan(BP_post1nan);
        GP = permute(GP_all(bands,:,:),[1,3,2]);
    case 'none'
        GP = true(nB,1,nCall);
    otherwise 
        error('optBP=%s is not defined.',optBP);
end
BP = isnan(GP);

%% read ADR transmission data
switch t_mode
    case {1,2,3}
        [ at_trans ] = load_ADR_VS('BINNING',WAdata.prop.binning,...
                                   'WAVELENGTH_FILTER',WAdata.prop.wavelength_filter);
    case {4}
        [ at_trans ] = load_T();
    otherwise
        error('Undefined t_mode %d',t_mode);
end
T = at_trans(:,:,bands); T(T<=1e-8) = nan;
T = permute(T,[3,1,2]); logT = log(T);

fprintf('finish loading ADR\n');
%%
% clear variables no longer used
clear Yif T at_trans BPdata1 BPdata2 BPdata_post BP_pri1nan BP_pri

%% main loop
if strcmpi(precision,'single')
    logYif = single(logYif);
    logT   = single(logT);
    WAb     = single(WAb);
end

% evaluate valid columns
if isempty(column_idxes), column_idxes = 1:nCall; end
valid_columns = find(~all(isnan(WAb),1));
Columns_valid = intersect(column_idxes,valid_columns);
if isempty(Columns_valid), fprintf(2,'Specified columns are invalid\n'); end

logT_extrap = logT;
for c = Columns_valid
    if any(isnan(logT(:,:,c)))
        if ~any(isnan(logT(:,:,c+1)))
            logT_extrap(:,:,c) = logT(:,:,c+1);
        elseif ~any(isnan(logT(:,:,c-1)))
            logT_extrap(:,:,c) = logT(:,:,c-1);
        elseif ~any(isnan(logT(:,:,c+2)))
            logT_extrap(:,:,c) = logT(:,:,c+2);
        elseif ~any(isnan(logT(:,:,c-2)))
            logT_extrap(:,:,c) = logT(:,:,c-2);
        else
            logT_extrap(:,:,c) = logT(:,:,c);
        end
    else
        logT_extrap(:,:,c) = logT(:,:,c);
    end
end

switch verbose
    case 0
        verbose_lad = 'no';
        verbose_huwacb = 'no';
        debug_huwacb = false;
        debug_lad = false;
    case 1
        verbose_lad = 'yes';
        verbose_huwacb = 'yes';
        debug_huwacb = false;
        debug_lad = false;
    case 2
        verbose_lad = 'yes';
        verbose_huwacb = 'yes';
        debug_huwacb = true;
        debug_lad = true;
    otherwise
        error('VIS=%d is not defined',verbose);
end

fprintf('Start processing\n');
tstart = datetime('now','TimeZone','America/New_York','Format','d-MMM-y HH:mm:ss Z');
fprintf('Current time is %s.\n',tstart);

switch upper(PROC_MODE)
    case {'CPU_1','GPU_1'}
        %%
        Yif_cor = nan([nLall,nCall,nBall],precision);
        Yif_isnan = nan([nLall,nCall,nBall],precision);
        T_est = nan([nBall,nCall],precision);
        Bg_est = nan([nLall,nCall,nBall],precision);
        AB_est = nan([nLall,nCall,nBall],precision);
        if ~isempty(opticelib)
            Ice_est = nan([nLall,nCall,nBall],precision);
        end
        ancillaries = struct('X',cell(1,nCall),'gp_bool',cell(1,nCall));
        Yif_cor_ori = nan([nLall,nCall,nBall],precision);
        Valid_pixels = false([nLall,nCall]);
        
        for c = Columns_valid
            tic;
            if Alib_out
                [Alib,infoAall,valid_idx] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,...
                bands_opt,WAb(:,c),cntRmvl);
            else
                [Alib] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,...
                bands_opt,WAb(:,c),cntRmvl);
            end
            if ~isempty(opticelib)
                if Alib_out
                    [Aicelib,infoAiceall,valid_idx_ice] = loadlibc_crism_icelib(opticelib,basenameWA,c,...
                        bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                else
                    [Aicelib] = loadlibc_crism_icelib(opticelib,basenameWA,c,...
                        bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                end
            else
                Aicelib = [];
            end
            if strcmpi(precision,'single')
                Alib = single(Alib); Aicelib = single(Aicelib);
            end

            [ logt_est,logYifc_cor,logAB,logBg,logIce,logYifc_cor_ori,logYifc_isnan,ancillary,vldpxl_c]...
                = sabcondc_v3l1_pub(Alib,logYif(:,:,c),WAb(:,c),logT_extrap(:,:,c),...
                  'GP',GP(:,:,c),...
                  'LAMBDA_A',lambda_a,'NITER',nIter,'PRECISION',precision,'GPU',gpu,...
                  'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                  'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                  'Aicelib',Aicelib);

            Yif_cor(lBool,c,bBool) = reshape(logYifc_cor',[nL,1,nB]);
            Yif_cor_ori(lBool,c,bBool) = reshape(logYifc_cor_ori',[nL,1,nB]);
            T_est(bBool,c) = logt_est;
            Bg_est(lBool,c,bBool) = reshape(logBg',[nL,1,nB]);
            AB_est(lBool,c,bBool) = reshape(logAB',[nL,1,nB]);
            Yif_isnan(lBool,c,bBool) = reshape(logYifc_isnan',[nL,1,nB]);
            if ~isempty(opticelib)
                Ice_est(lBool,c,bBool) = reshape(logIce',[nL,1,nB]);
            end
            ancillaries(c) = ancillary;
            Valid_pixels(lBool,c) = vldpxl_c';
            
            if Alib_out
                ancillaries(c).Alib = Alib;
                ancillaries(c).infoA = infoAall(valid_idx);
                if ~isempty(opticelib)
                    ancillaries(c).Aicelib = Aicelib;
                    ancillaries(c).infoAicelib = infoAiceall(valid_idx_ice);
                end
            end
            

            toc;
        end

        Yif_cor = exp(Yif_cor); T_est = exp(T_est);
        Bg_est = exp(Bg_est); AB_est = exp(AB_est);
        Yif_cor_ori = exp(Yif_cor_ori); 
        if ~isempty(opticelib)
            Ice_est = exp(Ice_est);
        end
    
    case {'CPU_2','GPU_2','GPU_BATCH_2'}
        %%
        Yif_cor = nan([nBall,nLall,nCall],precision);
        Yif_isnan = nan([nBall,nLall,nCall]);
        T_est = nan([nBall,1,nCall],precision);
        Bg_est = nan([nBall,nLall,nCall],precision);
        AB_est = nan([nBall,nLall,nCall],precision);
        badspcs = true([1,nLall,nCall]);
        Ice_est = nan([nBall,nLall,nCall],precision);
        Yif_cor_ori = nan([nBall,nLall,nCall],precision);
        X = [];
        
        n_batch = ceil(length(Columns_valid)/batch_size);
        
        for ni = 1:n_batch
            if ni~=n_batch
                Columns = Columns_valid((1+batch_size*(ni-1)):(batch_size*ni));
            elseif ni==n_batch
                Columns = Columns_valid((1+batch_size*(ni-1)):length(Columns_valid));
            end
            for i = 1:length(Columns)
                c = Columns(i);
                [Alib] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,...
                    bands_opt,WAb(:,c),cntRmvl);
                if ~isempty(opticelib)
                [Aicelib] = loadlibc_crism_icelib(opticelib,basenameWA,c,...
                    bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                else
                    Aicelib = [];
                end
                NA = size(Alib,2);
                if i==1
                    Alibs = Alib;
                    Aicelibs = Aicelib;
                else
                    Alibs = cat(3,Alibs,Alib);
                    Aicelibs = cat(3,Aicelibs,Aicelib);
                end
            end
            if strcmpi(precision,'single')
                Alibs = single(Alibs); Aicelibs = single(Aicelibs);
            end
            tic;
            [Yif_cor(bBool,lBool,Columns),T_est(bBool,1,Columns),...
                AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns),...
                Yif_isnan(bBool,lBool,Columns),Xc,badspcs(1,lBool,Columns)]...
            = sabcondc_v3l1_gpu_batch(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                      logT_extrap(:,:,Columns),...
                      BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                      'Aicelib',Aicelibs,'nIter',nIter,...
                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                      'gpu',gpu);
            switch upper(PROC_MODE)
                case {'GPU_BATCH_2'}
                Yif_cor_ori(bBool,lBool,Columns)...
                    = logYif(:,:,Columns) ...
                    - gather(pagefun(@mtimes,gpuArray(T_est(bBool,1,Columns)),gpuArray(Xc(1,:,:))))...
                    -Ice_est(bBool,lBool,Columns);
                otherwise
                    Yif_cor_ori(bBool,lBool,Columns) ...
                    = logYif(:,:,Columns)-T_est(bBool,1,Columns)*Xc(1,:,:)...
                      -Ice_est(bBool,lBool,Columns);
            end
            if ni==1
                X = nan(NA+1,nLall,nCall);
            end
            X(:,lBool,Columns) = Xc;
            toc;
        end
        Yif_cor = permute(Yif_cor,[2,3,1]);
        Yif_cor_ori = permute(Yif_cor_ori,[2,3,1]);
        Yif_isnan = permute(Yif_isnan,[2,3,1]);
        Bg_est = permute(Bg_est,[2,3,1]);
        AB_est = permute(AB_est,[2,3,1]);
        Ice_est = permute(Ice_est,[2,3,1]);
        T_est  = squeeze(T_est);
        X = permute(X,[2,3,1]);
        badspcs = squeeze(badspcs);
        Valid_pixels = ~badspcs;
    
        Yif_cor = exp(Yif_cor); T_est = exp(T_est);
        Bg_est = exp(Bg_est); AB_est = exp(AB_est); 
        Ice_est = exp(Ice_est);
        Yif_cor_ori = exp(Yif_cor_ori);
    otherwise
        error('Undefined PROC_MODE=%s',PROC_MODE);
        
end

tend = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
fprintf('finish procceing, current time is %s.\n',tend);
fprintf('Processing time is %s\n',tend-tstart);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Now saving...\n');

%% Write a setting file.
fname = [basename_cr '_settings.txt'];
if save_file
    fid = fopen(joinPath(save_dir,fname),'w');
else
    fid = 1; % standard output
end
hostname = char(java.net.InetAddress.getLocalHost.getHostName);
dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
fprintf(fid,'Hostname: %s\n',hostname);
fprintf(fid,'Username: %s\n',username);
fprintf(fid,'Time finished: %s\n',dt);

% ## I/O OPTIONS #---------------------------------------------------------
fprintf(fid,'SAVE_PDIR: %s\n',save_pdir);
fprintf(fid,'SAVE_DIR_YYYY_DOY: %d\n',save_dir_yyyy_doy);
fprintf(fid,'FORCE: %d\n',force);
fprintf(fid,'SKIP_IFEXIST: %d\n',skip_ifexist);
fprintf(fid,'ADDITIONAL_SUFFIX: %s\n',additional_suffix);
fprintf(fid,'INTERLEAVE_OUT: %s\n',interleave_out);
fprintf(fid,'SUBSET_COLUMNS_OUT: %d\n', subset_columns_out);
fprintf(fid,'ALIB_OUT: %d\n',Alib_out);

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
fprintf(fid,'OPT_IMG: %s\n',opt_img);
fprintf(fid,'TRRY_PDIR: %s\n',dir_yuk);
fprintf(fid,'FFC_IF_COUNTER: %d\n',ffc_counter);
fprintf(fid,'BANDS_OPT: %d\n',bands_opt);
fprintf(fid,'LINES:'); fprintf(fid,' %d',line_idxes); fprintf(fid, '\n');
fprintf(fid,'COLUMNS:'); fprintf(fid,' %d',column_idxes); fprintf(fid, '\n');
fprintf(fid,'METHODTYPE: %s\n',mt);
fprintf(fid,'OPTBP: %s\n',optBP);
fprintf(fid,'VERBOSE: %d\n', verbose);

% ## TRANSMISSION SPECTRUM OPTIONS #---------------------------------------
fprintf(fid,'T_MODE: %d\n', t_mode);

% ## LIBRARY OPTIONS #-----------------------------------------------------
fprintf(fid,'CNTRMVL: %d\n',cntRmvl);
fprintf(fid,'OPTINTERPID: %d\n',optInterpid);
fprintf(fid,'OPT_CRISMSPCLIB: %d\n',optCRISMspclib);
fprintf(fid,'OPT_RELAB: %d\n',optRELAB);
fprintf(fid,'OPT_SPLIBUSGS: %d\n',optUSGSsplib);
fprintf(fid,'OPT_CRISMTYPELIB: %d\n',optCRISMTypeLib);
fprintf(fid,'OPT_ICELIB: %d\n',opticelib);

% ## SABCONDC OPTIONS #----------------------------------------------------
fprintf(fid,'NITER: %d\n',nIter);
fprintf(fid,'LAMBDA_A:'); fprintf(fid,' %f',lambda_a); fprintf(fid,'\n');

% ## PROCESSING OPTIONS #--------------------------------------------------
fprintf(fid,'PRECISION: %s\n',precision);
fprintf(fid,'PROC_MODE: %s\n',PROC_MODE);

% ## ETCETERA #------------------------------------------------------------
% fprintf(fid,'GAUSSSIGMA: %f\n',gausssigma);

if fid>1
    fclose(fid);
end

%%
if save_file
    % hdr, mimics CAT file production
    hdr_cr = crism_const_cathdr(TRRIFdata,true,'DATE_TIME',dt);
    hdr_cr.cat_history = suffix;
    switch opt_img
        case 'if'
            hdr_cr.cat_input_files = [TRRIFdata.basename '.IMG'];
        case 'ra_if'
            hdr_cr.cat_input_files = [TRRRAIFdata.basename '_IF.IMG'];
        case {'TRRY','TRRC','TRRB'}
            hdr_cr.cat_input_files = [basenameTRRY '.IMG'];
        otherwise
            error('opt_img = %s is not defined',opt_img);
    end

    %% saving
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr '.hdr']));
    envihdrwritex(hdr_cr,joinPath(save_dir,[basename_cr '.hdr']),'OPT_CMOUT',false);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr '.img']));
    envidatawrite(single(Yif_cor),joinPath(save_dir,[basename_cr '.img']),hdr_cr);
    fprintf('Done\n');

    basename_ori = [basename_cr '_ori'];
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_ori '.hdr']));
    envihdrwritex(hdr_cr,joinPath(save_dir,[basename_ori '.hdr']),'OPT_CMOUT',false);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_ori '.img']));
    envidatawrite(single(Yif_cor_ori),joinPath(save_dir,[basename_ori '.img']),hdr_cr);
    fprintf('Done\n');

    fname_supple = joinPath(save_dir,[basename_cr '.mat']);
    wa = WAdata.img;
    wa = squeeze(wa)';
    fprintf('Saving %s ...\n',fname_supple);
    switch PROC_MODE
        case {'CPU_1','GPU_1'}
            save(fname_supple,'wa','bands','line_idxes','T_est','ancillaries','Valid_pixels');
        case {'GPU_BATCH_2'}
            save(fname_supple,'wa','bands','line_idxes','T_est','X','Yif_isnan');
    end
    fprintf('Done\n');

    basename_Bg = [basename_cr '_Bg'];
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg '.hdr']));
    envihdrwritex(hdr_cr,joinPath(save_dir,[basename_Bg '.hdr']),'OPT_CMOUT',false);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg '.img']));
    envidatawrite(single(Bg_est),joinPath(save_dir, [basename_Bg '.img']),hdr_cr);
    fprintf('Done\n');

    basename_AB = [basename_cr '_AB'];
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB '.hdr']));
    envihdrwritex(hdr_cr,joinPath(save_dir, [basename_AB '.hdr']),'OPT_CMOUT',false);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB '.img']));
    envidatawrite(single(AB_est),joinPath(save_dir, [basename_AB '.img']),hdr_cr);
    fprintf('Done\n');

    if ~isempty(opticelib)
        basename_Ice = [basename_cr '_Ice'];
        fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Ice '.hdr']));
        envihdrwritex(hdr_cr,joinPath(save_dir, [basename_Ice '.hdr']),'OPT_CMOUT',false);
        fprintf('Done\n');
        fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Ice '.img']));
        envidatawrite(single(Ice_est),joinPath(save_dir, [basename_Ice '.img']),hdr_cr);
        fprintf('Done\n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performing interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yif_cor_nr = Yif_cor;
Yif_mdl = AB_est .* Bg_est;
nan_cells = isnan(Yif_cor);
Yif_cor_nr(nan_cells) = Yif_mdl(nan_cells);

Valid_pixels_good = double(Valid_pixels);
Valid_pixels_good(Valid_pixels==0) = nan;
for bi=1:nBall
    Yif_cor_nr(:,:,bi) = Yif_cor_nr(:,:,bi) .* Valid_pixels_good;
end

% replace NaN with interpolation from a model
if save_file
    basename_cr_nr = [basename_cr '_nr'];
    hdr_cr_nr = hdr_cr;
    dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
    hdr_cr_nr.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
    hdr_cr_nr.cat_history = [hdr_cr_nr.cat_history '_nr'];
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr '.hdr']));
    envihdrwritex(hdr_cr_nr,joinPath(save_dir,[basename_cr_nr '.hdr']),'OPT_CMOUT',false);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr '.img']));
    envidatawrite(single(Yif_cor_nr),joinPath(save_dir,[basename_cr_nr '.img']),hdr_cr);
    fprintf('Done\n');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % gaussian filter
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logimg_sabcondl1_nan_replaced_ab = log(img_sabcondl1_nan_replaced) - log(Bg_est);
% logimg_gfl = nan(size(logimg_sabcondl1_nan_replaced_ab));
% gausswidth = gausssigma; fltsize = 5;
% h = fspecial('gaussian',fltsize,gausswidth);
% for i=1:size(img_sabcondl1_nan_replaced,3)
%     logimg_gfl(:,:,i) = nanimfilter(logimg_sabcondl1_nan_replaced_ab(:,:,i),h);
% end
% img_sabcondl1_nr_gf = exp(logimg_gfl) .* Bg_est;
% 
% basename_cr_nr_gf = [basename_cr_nr '_gf'];
% hdr_cr_nr_gf = hdr_cr_nr;
% dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
% hdr_cr_nr_gf.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced and gauss filtered after processing modified.}',dt);
% hdr_cr_nr_gf.cat_history = [hdr_cr_nr_gf.cat_history '_gf'];
% hdr_cr_nr_gf.gauss_filter_std = gausswidth;
% hdr_cr_nr_gf.gauss_filter_size = fltsize;
% hdr_cr_nr.cat_input_files = basename_cr_nr;
% 
% fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr_gf '.hdr']));
% envihdrwritex(hdr_cr_nr_gf,joinPath(save_dir,[basename_cr_nr_gf '.hdr']),'OPT_CMOUT',false);
% fprintf('Done\n');
% fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr_gf '.img']));
% envidatawrite(single(img_sabcondl1_nr_gf),joinPath(save_dir,[basename_cr_nr_gf '.img']),hdr_cr_nr_gf);
% fprintf('Done\n');

fprintf('Process completed!\n');

if save_file, diary off; end

%% Construct output
out = [];
if nargout==0
    
elseif nargout==1
    WA = squeeze(WAdata.img(:,:,:))';
    % take the subset of the columns
    if subset_columns_out
        Yif_cor     = Yif_cor(:,column_idxes,:);
        Yif_cor_nr  = Yif_cor_nr(:,column_idxes,:);
        Yif_cor_ori = Yif_cor_ori(:,column_idxes,:);
        Yif_isnan   = Yif_isnan(:,column_idxes,:);
        AB_est      = AB_est(:,column_idxes,:);
        Bg_est      = Bg_est(:,column_idxes,:);
        WA          = WA(:,column_idxes);
        T_est       = T_est(:,column_idxes);
        Valid_pixels= Valid_pixels(:,column_idxes);
        if ~isempty(opticelib)
            Ice_est = Ice_est(:,column_idxes,:);
        end
    end
    % Permute the output.
    prmt_ordr   = [find(interleave_out(1)==interleave_default),...
                   find(interleave_out(2)==interleave_default),...
                   find(interleave_out(3)==interleave_default)];
    Yif_cor     = permute(Yif_cor,    prmt_ordr);
    Yif_cor_nr  = permute(Yif_cor_nr, prmt_ordr);
    Yif_cor_ori = permute(Yif_cor_ori,prmt_ordr);
    Yif_isnan   = permute(Yif_isnan,  prmt_ordr);
    AB_est      = permute(AB_est,     prmt_ordr);
    Bg_est      = permute(Bg_est,     prmt_ordr);
    Valid_pixels= permute(Valid_pixels,prmt_ordr);
    if ~isempty(opticelib)
        Ice_est = permute(Ice_est,    prmt_ordr);
    end

    out.Yif_cor     = Yif_cor;
    out.Yif_cor_nr  = Yif_cor_nr;
    out.Yif_cor_ori = Yif_cor_ori;
    out.AB_est      = AB_est;
    out.Bg_est      = Bg_est;
    if ~isempty(opticelib)
        out.Ice_est = Ice_est;
    end
    out.T_est        = T_est;
    out.Yif_isnan    = Yif_isnan;
    out.Valid_pixels = Valid_pixels;
    out.WA           = WA;
    out.lines        = line_idxes;
    out.columns      = column_idxes;
    out.bands        = bands;
    out.interleave_out     = interleave_out;
    out.subset_columns_out = subset_columns_out;
    
    switch upper(PROC_MODE)
        case {'CPU_1','GPU_1'}
            if subset_columns_out
                ancillaries  = ancillaries(column_idxes);
            end
            out.ancillaries  = ancillaries;
        case {'CPU_2','GPU_2','GPU_BATCH_2'}
            X = permute(X(:,column_idxes,:),prmt_ordr);
            out.X = X;
        otherwise
            error('Undefined PROC_MODE=%s',PROC_MODE);
    end
end


end
