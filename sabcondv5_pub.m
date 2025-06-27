function [out] = sabcondv5_pub(obs_id,varargin)
% [out] = sabcondv5_pub(obs_id,varargin)
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
%           out.obs_id      = obs_id;
%           out.TRR3IF      = TRRIFdata.basename;
%           out.Yif_cor     = Yif_cor;
%           out.Yif_cor_nr  = Yif_cor_nr;
%           out.Yif_cor_ori = Yif_cor_ori;
%           out.Yif_nr_ds   = Yif_nr_ds;
%           out.AB_est      = AB_est;
%           out.Bg_est      = Bg_est;
%           out.AB_est_ds   = AB_est_ds;
%           out.Bg_est_ds   = Bg_est_ds;
%           out.Yifmdl_ds   = Yifmdl_ds;
%           out.T_est       = T_est;
%           out.Yif_isnan    = Yif_isnan;
%           out.Valid_pixels = Valid_pixels;
%           out.WA          = WA;
%           out.WV          = WV;
%           out.lines       = line_idxes;
%           out.columns     = Columns;
%           out.bands       = bands;
%           out.interleave_out     = interleave_out;
%           out.subset_columns_out = subset_columns_out;
%           out.GP = GP;
%           out.BP = BP;
%           out.ancillaries = ancillaries;
%
% OPTIONAL Parameters 
%  ## I/O OPTIONS  #-------------------------------------------------------
%   'SAVE_FILE': boolean
%       whether or not to save processed images. If true, two optioal 
%       parameters 'FORCE','SKIP_IFEXIST' have no effect.
%       (default) true
%   'STORAGE_SAVING_LEVEL': string
%       determine how much to save storage, no effect under save_file=1
%       Normal  - All the byproducts are saved
%       Highest - Only nr_ds and mdl_ds are saved, bands are also scropped.
%       (default) Normal
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
%       (default) ''
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
%   'CROP_BANDS': Boolean
%       whether or not to crop unprocessed bands that are filled with NaN
%       by default. Set to true if you set STORAGE_SAVING_LEVEL=Highest
%       (default) false
%
%  ## INPUT IMAGE OPTIONS #------------------------------------------------
%   'OPT_IMG': string, {'IF','RA_IF','TRRY','TRRB','TRRC'}
%       type of input image to be used
%       (default) 'TRRB'
%   'IMG_CUBE': image cube [L x S x B]
%       if not empty, this input image will be used. if empty, Image 
%       corresponding to 'OPT_IMG' is read
%       from file. 
%       (default) []
%   'IMG_CUBE_BAND_INVERSE': boolean
%       need this when you input 'IMG_CUBE' to clarify the image.
%       (default) []
%   'TRRY_PDIR': any string
%       root directory path where {'TRRY','TRRB','TRRC'} images are stored.
%       If you choose either {'IF', 'RA_IF'} for 'OPT_IMG', you do not need
%       to specify this.
%   'FFC_IF_COUNTER': integer, {1,3}
%       Observation counter of the image to be processed. Used only for 
%       processing FFC images. 
%       (default) 1
%   'OBS_COUNTER_SCENE': string
%        regular expression to match observation counter, such as
%        '0[13]{1}', observation counter for scene measurements.
%        (default) depends on observation type
%   'OBS_COUNTER_DF': 
%        regular expression to match observation counter, such as
%        '0[13]{1}', observation counter for dark frame measurements.
%        (default) depends on observation type
%
%  ## GENERAL SABCOND OPTIONS #--------------------------------------------
%   'BANDS_OPT' : integer, {4}
%       an option for wavelength channels to use. This is the 
%       input for crmsab_genBands()
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
%   'COLUMN_SKIP: integer
%       the stride for which processing columns will be skipped.
%       (default) 1
%   'WEIGHT_MODE': integer {0,1}
%       mode id of weight option. 0 for uniform weight, 1 for weight based
%       on noise estimation.
%       (default) 0
%   'LAMBDA_UPDATE_RULE': string {'L1SUM','MED','NONE'}
%       define how to update trade-off parameters.
%       (default) 'L1SUM'
%   'THRESHOLD_BADSPC': scalar
%       threshold value for which the spectra is considered to be
%       completely corrupted
%       (default) 0.8
%   'FFC_MODE': boolean
%       whether or not to perform FFC correcion
%       (default) false
%   'OPT_BANDS_IGNORE_INIT': string, {'none','ltn035'}
%       options for selecting bands to be ignored in the first iteration
%       (default) 'none'
%   'INCLUDE_ICE': bool,
%       whether or not to include the contribution of ICE from (Aicelib)to 
%       the model
%       (default) false
%
%  ## PRE-PROCESSING OPTIONS #---------------------------------------------
%   'CAL_BIAS_COR': Integer {0,1,2}
%       Option for how to perform image based bias correction.
%       0: none, 1: conservative correction, 2: aggressive correction
%       (default) 0
%
%  ## TRANSMISSION SPECTRUM OPTIONS #--------------------------------------
%    'T_MODE': integer {1,2,3,4,5}
%       option for what kind of transmission spectrum frame is used.
%       (default) 2
%    'OBS_ID_T': string
%       obs_id of the transmission file used for processing. Only valid for
%       T_MODE=5.
%       (default) ''
%    'VARARGIN_T': cell
%       varargin for the function 'crmsab_load_T_given' or 'crmsab_load_T_sclk_closest'
%       (default) {}
%
%  ## LIBRARY OPTIONS #----------------------------------------------------
%   'OPT_LIBRARY': char, string {'mini','full'}
%       'mini': if you have the minimum set of the library for default
%       processing, you need to use band_opt= 4 and cntRmvl option is just 
%       ignored. library is also limited to 1112 and 1164
%       'full': full support. This is for custom setting for the library.
%       (default) 'mini'
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
%   'T_UPDATE: scalar,
%       how many times t is updated
%       (default) inf (update as many as nIter)
%   'LOGT_NEG': boolean,
%       whether or not to force logT to be negative
%       (default) false
%   'LOGT_RELAX': boolean,
%       whether or not to relax logT update with the residuals
%       (default) false
%
%  ## PROCESSING OPTIONS #-------------------------------------------------
%   'PRECISION': string, {'single','double'}
%       percision with which processing is performed.
%       (default) 'double'
%   'PROC_MODE': string, {'CPU_1','GPU_1','CPU_2','GPU_2','GPU_BATCH_2','CPU_3','GPU_3','GPU_BATCH_3'}
%       option for which processor is used CPU or GPU, and which algorithm
%       is used, formatted as PPP[_BATCH]_A, where PPP represents 
%       processor, and A is the index for algorithm. BATCH indicates 
%       (default) 'CPU_1'
%       **ALGORITHM description**
%           1: Bad pixel detected during processing are substitueded by
%              model values.
%           2: Bad pixel detected during processing are exactly ignored
%              during the processing.
%           3: ONLY multiplicative bias constant is learned instead of
%              transmission spectrum
%   'BATCH_SIZE': integer
%       column size for which processing is performed. Valid only if
%       'GPU_BATCH_*' mode is selected.
%       (default) 10
%   'DEBUG': boolean
%        going into debug mode
%        (default) false

global crism_env_vars

% ## I/O OPTIONS #---------------------------------------------------------
save_file          = true;
storage_saving_level = 'NORMAL';
save_pdir          = './resu/';
save_dir_yyyy_doy  = false;
force              = false;
skip_ifexist       = false;
additional_suffix  = '';
interleave_out     = 'lsb';
interleave_default = 'lsb';
subset_columns_out = false;
Alib_out           = false;
do_crop_bands      = false;

% ## INPUT IMAGE OPTIONS #-------------------------------------------------
opt_img      = 'TRRB';
img_cube     = [];
img_cube_band_inverse = [];
dir_TRRX      = crism_env_vars.dir_TRRX; % TRRY_PDIR
ffc_counter  = 1;
OBS_COUNTER_SCENE_custom = 0;
OBS_COUNTER_DF_custom = 0;
take_median = false;

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
bands_opt    = 4;
line_idxes   = [];                     % LINES
column_idxes = [];
mt           = 'sabcondpub_v1';        % METHODTYPE
optBP        = 'pri';                  %{'pri','all','none'}
verbose      = 0;
column_skip  = 1;
weight_mode  = 0;
lambda_update_rule = 'L1SUM';
th_badspc    = 0.5;
ffc_mode     = false;
opt_bands_ignore_init = 'none';
include_ice = false;

% ## PRE-PROCESSING OPTIONS #----------------------------------------------
cal_bias_cor = 0;

% ## TRANSMISSION SPECTRUM OPTIONS #---------------------------------------
t_mode = 2;
obs_id_T = '';
varargin_T = {};

% ## LIBRARY OPTIONS #-----------------------------------------------------
optLibrary      = 'mini'; %{'mini','full'}
cntRmvl         = 1;
optInterpid     = 1;
optCRISMspclib  = 1;
optRELAB        = 1;
optUSGSsplib    = 1;
optCRISMTypeLib = 2;
opticelib       = '';
optsurficelib   = '';

% ## SABCONDC OPTIONS #----------------------------------------------------
nIter = 5;
lambda_a = 0.01;
lambda_a_ice   = 0;
lambda_a_surfice   = 0;
t_update     = inf;
logT_neg = false;
logt_relax = false;
% thrd_noise_small = 0.0015;
thrd_snr = 0.02;

% ## PROCESSING OPTIONS #--------------------------------------------------
precision  = 'double';
PROC_MODE  = 'CPU_1';
batch_size = 10;
is_debug   = false;
debug_l_plot = 1;

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
            case 'STORAGE_SAVING_LEVEL'
                storage_saving_level = varargin{i+1};
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
            case 'CROP_BANDS'
                do_crop_bands = varargin{i+1};
                
            % ## INPUT IMAGE OPTiONS #-------------------------------------
            case 'OPT_IMG'
                opt_img = varargin{i+1};
            case 'IMG_CUBE'
                img_cube = varargin{i+1};
            case 'IMG_CUBE_BAND_INVERSE'
                img_cube_band_inverse = varargin{i+1};
            case 'TRRY_PDIR'
                dir_TRRX = varargin{i+1};
            case 'FFC_IF_COUNTER'
                ffc_counter = varargin{i+1};
            case 'OBS_COUNTER_SCENE'
                obs_counter_tmp = varargin{i+1};
                OBS_COUNTER_SCENE_custom = 1;
            case 'OBS_COUNTER_DF'
                obs_counter_df_tmp = varargin{i+1};
                OBS_COUNTER_DF_custom = 1;
            case 'TAKE_MEDIAN'
                take_median = varargin{i+1};
                
            % ## GENERAL SABCOND OPTIONS #---------------------------------
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
            case 'COLUMN_SKIP'
                column_skip = varargin{i+1};
            case 'WEIGHT_MODE'
                weight_mode = varargin{i+1};
            case 'LAMBDA_UPDATE_RULE'
                lambda_update_rule = varargin{i+1};
            case 'THRESHOLD_BADSPC'
                th_badspc = varargin{i+1};
            case 'FFC_MODE'
                ffc_mode = varargin{i+1};
            case 'OPT_BANDS_IGNORE_INIT'
                opt_bands_ignore_init = varargin{i+1};
            case 'INCLUDE_ICE'
                include_ice = varargin{i+1};
                
            % % ## PRE-PROCESSING OPTIONS #--------------------------------
            case 'CAL_BIAS_COR'
                cal_bias_cor = varargin{i+1};
                
            % ## TRANSMISSION SPECTRUM OPTIONS #---------------------------
            case 'T_MODE'
                t_mode = varargin{i+1};
            case 'OBS_ID_T'
                obs_id_T = varargin{i+1};
            case 'VARARGIN_T'
                varargin_T = varargin{i+1};
                
            % ## LIBRARY OPTIONS #-----------------------------------------
            case 'LIBRARY_OPT'
                optLibrary = varargin{i+1};
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
            case 'OPT_SURFICELIB'
                optsurficelib = varargin{i+1};
                
            % ## SABCONDC OPTIONS #----------------------------------------
            case 'NITER'
                nIter = varargin{i+1};
            case 'LAMBDA_A'
                lambda_a = varargin{i+1};
            case 'LAMBDA_A_ICE'
                lambda_a_ice = varargin{i+1};
            case 'LAMBDA_A_SURFICE'
                lambda_a_surfice = varargin{i+1};
            case 'T_UPDATE'
                t_update = varargin{i+1};
            case 'LOGT_NEG'
                logT_neg = varargin{i+1};
            case 'LOGT_RELAX'
                logt_relax = varargin{i+1};
            % case 'THRESHOLD_NOISE_SMALL'
            %     thrd_noise_small = varargin{i+1};
            case 'THRESHOLD_SNR'
                thrd_snr = varargin{i+1};
                
            % ## PROCESSING OPTIONS #--------------------------------------
            case 'PRECISION'
                precision = varargin{i+1};
            case 'PROC_MODE'
                PROC_MODE = varargin{i+1};
            case 'BATCH_SIZE'
                batch_size = varargin{i+1};
            case 'DEBUG'
                is_debug = varargin{i+1};
            case 'DEBUG_L_PLOT'
                debug_l_plot = varargin{i+1};
                
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if save_file && force && skip_ifexist
    error('You are forcing or skipping? Not sure what you want');
end

if save_file && ~exist(save_pdir,'dir')
    [status] = mkdir(save_pdir); 
    if status
        if verbose, fprintf('"%s" is created.\n',save_pdir); end
        chmod777(save_pdir,verbose);
    else
        error('Failed to create %s',save_pdir);
    end
end

switch upper(storage_saving_level)
    case 'HIGHEST'
        if ~do_crop_bands
            if verbose
                fprintf('crop_bands is always set to true with STORAGE_SAVING_LEVEL=HIGHEST\n');
            end
        end
        do_crop_bands = true;
end


optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
% libprefix = crmsab_const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,'');

switch upper(PROC_MODE)
    case {'CPU_1','CPU_2','CPU_3','CPU_4'}
        gpu = false;
    case {'GPU_1','GPU_2','GPU_BATCH_2','GPU_3','GPU_BATCH_3','GPU_4','GPU_BATCH_4'}
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

% determine batch size for each PROC_MODE----------------------------------
switch upper(PROC_MODE)
    case {'CPU_1','CPU_2','GPU_2','CPU_3','GPU_3','CPU_4','GPU_4'}
        batch_size = 1;
    case {'GPU_BATCH_2','GPU_BATCH_3','GPU_BATCH_4'}
        
    otherwise
        error('Undefined PROC_MODE=%s',PROC_MODE);
end

%% INPUT Image options
[ yyyy_doy,obs_classType ] = crism_searchOBSID2YYYY_DOY_v2(obs_id);
switch obs_classType
    case {'FRT','HRL','HRS'}
        obs_counter = '07';
        obs_counter_epf = '[0-689A-Za-z]{2}';
        obs_counter_epfdf = '0[0E]{1}';
        obs_counter_df = '0[68]{1}';
    case {'FRS','ATO'}
        obs_counter = '01';
        obs_counter_df = '0[03]{1}';
        obs_counter_epf = '';
        obs_counter_epfdf = '';
        obs_counter_un = '02';
    case 'FFC'
        obs_counter = '0[13]{1}';
        obs_counter_df = '0[024]{1}';
        % this could be switched.
        obs_counter_epf = '';
    case 'CAL'
        obs_counter = '[0-9a-fA-F]{2}';
        obs_counter_df = '[0-9a-fA-F]{1}';
    case 'ICL'
        obs_counter = '[0-9a-fA-F]{2}';
        obs_counter_df = '[0-9a-fA-F]{1}';
    case {'MSP','HSP'}
        obs_counter = '01';
        obs_counter_df = '0[02]{1}';
        obs_counter_epf = '';
        obs_counter_epfdf = '';
    otherwise
        error('OBS_TYPE %s is not supported yet.',obs_classType);
end

if OBS_COUNTER_SCENE_custom
    obs_counter = obs_counter_tmp;
end
if OBS_COUNTER_DF_custom
    obs_counter_df = obs_counter_df_tmp;
end

%% Read image and ancillary data and format them for processing
% crism_obs_info = crism_get_obs_info_v2(obs_id, 'SENSOR_ID', 'L');
crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L',...
    'obs_counter_scene',obs_counter,'obs_counter_df',obs_counter_df);
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

[DFdata1,DFdata2] = crism_get_DFdata4SC(TRRIFdata,crism_obs);

%%
%-------------------------------------------------------------------------%
% cheking the file exist or not.
%-------------------------------------------------------------------------%
if save_file
    if save_dir_yyyy_doy
        dirpath_yyyy_doy = joinPath(save_pdir,crism_obs.info.yyyy_doy);
        if ~exist(dirpath_yyyy_doy,'dir')
            status = mkdir(dirpath_yyyy_doy);
            if status
                if verbose, fprintf('"%s" is created.\n',dirpath_yyyy_doy); end
                chmod777(dirpath_yyyy_doy,verbose);
            else
                error('Failed to create %s',dirpath_yyyy_doy);
            end
        end
        save_dir = joinPath(dirpath_yyyy_doy,crism_obs.info.dirname);
        if ~exist(save_dir,'dir')
            status = mkdir(save_dir);
            if status
                if verbose, fprintf('"%s" is created.\n',save_dir); end
                chmod777(save_dir,verbose);
            else
                error('Failed to create %s',save_dir);
            end
        end
    else
        save_dir = joinPath(save_pdir,crism_obs.info.dirname);
        if ~exist(save_dir,'dir')
            status = mkdir(save_dir);
            if status
                if verbose, fprintf('"%s" is created.\n',save_dir); end
                chmod777(save_dir,verbose);
            else
                error('Failed to create %s',save_dir);
            end
        end
    end
end

suffix = crmsab_const_suffix_v3(mt,additional_suffix);

fprintf('suffix will be \n"%s"\n',suffix);

switch upper(opt_img)
    case 'IF'
        basename_cr = [TRRIFdata.basenameIF suffix];
    case 'RA_IF'
        basename_cr = [crism_obs.info.basenameRA '_IF' suffix];
    case {'TRRY','TRRB','TRRC','TRRD','TRRE'}
        trr_vr = opt_img(4);
        basenameTRRY = crism_get_TRRXbasename(TRRIFdata,trr_vr);
        basename_cr = [basenameTRRY suffix];
    otherwise
        error('opt_img = %s is not defined',opt_img);
end


if save_file 
    fpath_cr = joinPath(save_dir,[basename_cr,'_nr_ds.img']);
    if exist(fpath_cr,'file')
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
end


if save_file && ~exist(save_dir,'dir'), mkdir(save_dir); end

% open log file
username = char(java.lang.System.getProperty('user.name'));

if save_file
    fname_log = sprintf('log_%s_%s.txt',username,datetime('now','TimeZone','local','Format','yyyyMMdd'));
    fpath_diary = joinPath(save_dir,fname_log);
    diary(fpath_diary);
end

%% Read image and ancillary data and format them for processing
TRRIFdata.load_basenamesCDR();
WAdata = TRRIFdata.readCDR('WA'); WAdata.readimgi();
bands = crmsab_genBands_v2(WAdata.prop.wavelength_filter,bands_opt,WAdata.prop.binning,WAdata.prop.sclk);


nLall = TRRIFdata.hdr.lines; nCall = TRRIFdata.hdr.samples; nBall = TRRIFdata.hdr.bands;
if isempty(line_idxes), line_idxes = 1:nLall; end
nL = length(line_idxes); nB = length(bands);
lBool = false(nLall,1); lBool(line_idxes) = true;
bBool = false(nBall,1); bBool(bands) = true;

basenameWA = WAdata.basename;
WAb = squeeze(WAdata.img(:,:,bands))';

if isempty(img_cube)
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

        case {'TRRY','TRRB','TRRC','TRRD','TRRE'}
            d_IoF = joinPath(dir_TRRX, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
            TRRYIFdata = CRISMdata(basenameTRRY,d_IoF);
            Yif = TRRYIFdata.readimgi();

        otherwise
            error('opt_img = %s is not defined',opt_img);
    end
    img_cube_isempty = true;
else
    % You will perform log conversion, so cast them to double precision.
    Yif = double(img_cube);
    % Check the size of the image just in case.
    [Lcube,Scube,Bcube] = size(Yif);
    if Lcube~=nLall || Scube~=nCall || Bcube~=nBall
        error('The size of "IMG_CUBE" is wrong');
    end
    % You need to specify the direction of the band of the image
    if isempty(img_cube_band_inverse)
        error('Specify "IMG_CUBE_BAND_INVERSE" when you enter "IMG_CUBE"');
    else
        if img_cube_band_inverse==0
            Yif = flip(Yif,3);
        end
    end
    img_cube_isempty = false;
    clear img_cube;
end


TRRIFdata.readCDR('DM');
DMmask = TRRIFdata.cdr.DM.readimgi();
fprintf('finish loading Image\n');

%% Compute calibration biases
% put nan for invalid values and extract lines focused.
Yif(Yif<=1e-8) = nan;
Yif = Yif(line_idxes,:,:);

% estimate potential biases
bands4bias = 1:252;
switch cal_bias_cor
    case {1,2}
        mat_name = [TRRYIFdata.basename sprintf('_BPpost3%dt%d.mat',bands4bias(1),bands4bias(end))];
        if exist(mat_name,'file')
            load(mat_name,'dev_coef','band_bias_std','bands4bias','dev_sub');
        else
            [BPpost_plus,band_bias_std,dev_sub,val_ratio_3d2,val_ratio_3d2_smth2] = detect_BPpost_wBias3(...
                Yif(:,:,:), DMmask,'bands',bands4bias,'debug',is_debug);
            
            [dev_coef,img_nanmed_estimate,dev_coef_ori] = calc_deviation_BPpost2(...
                BPpost_plus(1,:,bands4bias),Yif(:,:,bands4bias));
            save(mat_name,'dev_coef','band_bias_std','dev_sub','img_nanmed_estimate','dev_coef_ori','bands4bias','BPpost_plus','val_ratio_3d2','val_ratio_3d2_smth2');
        end
        % bands_bias_mad is for later processing.
        bands_bias_mad = permute(band_bias_std(:,:,bands),[3,1,2]) .* norminv(0.75);
        
        switch cal_bias_cor
            case 1 % conservative correction
                % dev_coef = dev_coef;
            case 2 % aggressive correction
                dev_coef = Yif(:,:,bands4bias) ./ (Yif(:,:,bands4bias)-dev_sub);
        end
        % apply biases for severely corrupted ones.
        Yif = Yif(:,:,bands4bias)./dev_coef;
        
    case {0}
        dev_coef = ones(1,size(Yif,2),length(bands4bias));
        bands_bias_mad = zeros(length(bands),1);
        % apply biases for severely corrupted ones.
        % Yif = Yif(:,:,bands4bias)./dev_coef;
    otherwise
        error('Undefined CAL_BIAS_COR=%d',cal_bias_cor);
end

Yif(516,:,1:48) = nan;
Yif(516,341:end,49) = nan;
Yif = Yif(:,:,bands);
Yif(Yif<=0) = nan;
logYif = log(Yif);
logYif = permute(logYif,[3,1,2]);

if take_median
    logYif = median(logYif, 2, 'omitnan');
    nLall = 1;
    lBool = true;
end

%%
% read bad pixel
if isempty(DFdata2)
    [BPdata1,BPdata2,BPdata_post] = crism_load_BPdataSC_fromDF(...
        TRRIFdata,DFdata1.basename,[]);
else
    [BPdata1,BPdata2,BPdata_post] = crism_load_BPdataSC_fromDF(...
        TRRIFdata,DFdata1.basename,DFdata2.basename);
end
switch lower(optBP)
    case 'pri'
        [BP1nanfull] = crism_formatBPpri1nan(BPdata1,BPdata2, ...
                            'band_inverse',true,'interleave','lsb');
        [GP1nanfull] = crism_convertBP1nan2GP1nan(BP1nanfull);
        GP1nan = permute(GP1nanfull(:,:,bands),[3,1,2]);
        BP1nan = permute(BP1nanfull(:,:,bands),[3,1,2]);
    case 'all'
        [BP1nanfull] = crism_formatBP1nan(BPdata_post, ...
                            'band_inverse',true,'interleave','lsb');
        [GP1nanfull] = crism_convertBP1nan2GP1nan(BP1nanfull);
        GP1nan = permute(GP1nanfull(:,:,bands),[3,1,2]);
        BP1nan = permute(BP1nanfull(:,:,bands),[3,1,2]);
    case 'none'
        GP1nan = ones(nB,1,nCall);
        BP1nan = nan(nB,1,nCall);
    otherwise 
        error('optBP=%s is not defined.',optBP);
end

BP = (BP1nan==1); GP = (GP1nan==1);

%% Weight mode
switch weight_mode
    case 0
        stdl1_ifdf  = zeros(nB,1,nCall,0);
        WA_um_pitch = zeros(nB,1,nCall,0);
        SFimg       = zeros(nB,1,nCall,0);
    case 1
        %% compute the weight for each dimension
        switch upper(opt_img)
            case {'IF','RA_IF'}
                error('weight_mode=%d only works with TRRB,TRRC,TRRC',weight_mode);
            case {'TRRY','TRRB','TRRC','TRRD','TRRE'}
                switch TRRIFdata.lbl.OBSERVATION_TYPE
                    case {'FRT','HRL','HRS','FFC'}
                        Noutliers = 4;
                    case {'FRS','ATO'}
                        Noutliers = 2;
                    otherwise
                        error('Please define for other cases')
                end
                % load processed dark files
                propDF1_IF = DFdata1.prop;
                propDF1_IF.activity_id = 'IF';
                propDF1_IF.product_type = 'TRR';
                propDF1_IF.version = trr_vr;
                bnameDF1_IF = crism_get_basenameOBS_fromProp(propDF1_IF);
                load(joinPath(d_IoF,[bnameDF1_IF '.mat']),'IoF_bk1_o');
                IoF_bk1_o = flip(IoF_bk1_o,3);
                IoF_bk1_o = permute(IoF_bk1_o,[3,1,2]);
                switch TRRIFdata.lbl.OBSERVATION_TYPE
                    case {'FRT','HRL','HRS','FFC'}
                        propDF2_IF = DFdata2.prop;
                        propDF2_IF.activity_id = 'IF';
                        propDF2_IF.product_type = 'TRR';
                        propDF2_IF.version = trr_vr;
                        bnameDF2_IF = crism_get_basenameOBS_fromProp(propDF2_IF);
                        load(joinPath(d_IoF,[bnameDF2_IF '.mat']),'IoF_bk2_o');
                        IoF_bk2_o = flip(IoF_bk2_o,3);
                        IoF_bk2_o = permute(IoF_bk2_o,[3,1,2]);
                    case {'FRS','ATO'}
                        IoF_bk2_o = [];
                    otherwise
                        error('Please define for other cases')
                end

                stdl1_ifdf = nan(size(GP));
                ifdfImage = cat(2,IoF_bk1_o,IoF_bk2_o);
                for cd = 1:nCall
                    ifdfimagec = ifdfImage(bands,:,cd);
                    % std_df = robust_v3('stdl1',ifdfimagec,2,'NOutliers',Noutliers);
                    std_df = robust_v3('med_abs_dev_from_med',ifdfimagec,2,'NOutliers',Noutliers);
                    stdl1_ifdf(:,:,cd) = std_df(:);
                end

                % read SFdata
                switch upper(opt_img)
                    case {'TRRY','TRRB'}
                        SFdata = TRRIFdata.readCDR('SF');

                    case {'TRRC','TRRD','TRRE'}
                        SFdata = TRRIFdata.readCDR('SF');

                    otherwise
                        error('opt_img = %s is not defined',opt_img);
                end

                % Estimate the amount of photon noise
                [WA_um_pitch] = crism_get_WA_um_pitch(WAdata);
                WA_um_pitch = permute(WA_um_pitch(:,:,bands),[3,1,2]);
                SFimg = SFdata.readimgi();
                SFimg = permute(SFimg(:,:,bands),[3,1,2]);
        end
    otherwise
        error('Undefined weight mode: %d',weight_mode);
end

%% read ADR transmission data
switch t_mode
    case {1,2,3}
        [ at_trans ] = crism_load_ADR_VS('BINNING',WAdata.prop.binning,...
            'WAVELENGTH_FILTER',WAdata.prop.wavelength_filter,'SCLK',WAdata.prop.sclk);
    case {4}
        sclk_img = (TRRIFdata.get_sclk_start()+TRRIFdata.get_sclk_stop())/2;
        [ at_trans ] = crmsab_load_T_sclk_closest(sclk_img,varargin_T{:});
        at_trans = crism_bin_image_frames(at_trans,'binning',WAdata.prop.binning);
    case {5}
        [ at_trans ] = crmsab_load_T_given_v2( obs_id_T,varargin_T{:});
        at_trans = crism_bin_image_frames(at_trans,'binning',WAdata.prop.binning);
    otherwise
        error('Undefined t_mode %d',t_mode);
end
T = at_trans(:,:,bands); T(T<=1e-8) = nan;
T = permute(T,[3,1,2]); 

switch precision
    case 'double'
        T = double(T);
end
logT = log(T);

fprintf('finish loading ADR\n');
%%
% clear variables no longer used
clear Yif T at_trans BPdata1 BPdata2 BPdata_post

%% main loop
if strcmpi(precision,'single')
    logYif = single(logYif);
    logT   = single(logT);
    WAb     = single(WAb);
end

% evaluate valid columns
if isempty(column_idxes), column_idxes = 1:nCall; end
column_idxes = column_idxes(1:column_skip:length(column_idxes));
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
        ancillaries = struct('Xt',cell(nCall,1),'Xlib',cell(nCall,1),...
                             'Xice',cell(nCall,1),...
                             'Alib',cell(nCall,1),'Aicelib',cell(nCall,1),...
                             'infoAlib',cell(nCall,1),'infoAicelib',cell(nCall,1));
        Yif_cor_ori = nan([nLall,nCall,nBall],precision);
        Valid_pixels = false([nLall,nCall]);
        
        for c = Columns_valid
            % tic;
            if Alib_out
                switch upper(optLibrary)
                    case 'FULL'
                        [Alib,infoAall,valid_idx] = crmsab_loadAlibconv(optLibs,basenameWA,optInterpid,c,...
                            bands_opt,WAb(:,c),cntRmvl);
                    case 'MINI'
                        [Alib,infoAall,valid_idx] = crmsab_loadAlibconv_mini(optLibs,basenameWA,optInterpid,c,bands_opt);
                    otherwise
                        error('Undefined OPT_LIBRARY %s',optLibrary);
                end
                ancillaries(c).Alib = Alib;
                ancillaries(c).infoAlib = infoAall(valid_idx);
            else
                switch upper(optLibrary)
                    case 'FULL'
                        [Alib] = crmsab_loadAlibconv(optLibs,basenameWA,optInterpid,c,...
                            bands_opt,WAb(:,c),cntRmvl);
                    case 'MINI'
                        [Alib] = crmsab_loadAlibconv_mini(optLibs,basenameWA,optInterpid,c,bands_opt);
                    otherwise
                        error('Undefined OPT_LIBRARY %s',optLibrary);
                end
            end
            if ~isempty(opticelib)
                if Alib_out
                    switch upper(optLibrary)
                        case 'FULL'
                            [Aicelib,infoAiceall,valid_idx_ice] = crmsab_loadAicelibconv(opticelib,basenameWA,c,...
                                bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                        case 'MINI'
                            [Aicelib,infoAiceall,valid_idx_ice] = crmsab_loadAicelibconv_mini(opticelib,basenameWA,c,bands_opt);
                        otherwise
                            error('Undefined OPT_LIBRARY %s',optLibrary);
                    end
                    ancillaries(c).Aicelib = Aicelib;
                    ancillaries(c).infoAicelib = infoAiceall(valid_idx_ice);
                else
                    switch upper(optLibrary)
                        case 'FULL'
                            [Aicelib] = crmsab_loadAicelibconv(opticelib,basenameWA,c,...
                                bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                        case 'MINI'
                            [Aicelib] = crmsab_loadAicelibconv_mini(opticelib,basenameWA,c,bands_opt);
                        otherwise
                            error('Undefined OPT_LIBRARY %s',optLibrary);
                    end
                end
            else
                Aicelib = [];
            end
            if ~isempty(optsurficelib)
                if Alib_out
                    switch upper(optLibrary)
                        case 'FULL'
                            [Asurficelib,infoAsurficeall,valid_idx_surfice] = crmsab_loadAsurficelibconv(optsurficelib,basenameWA,c,...
                                bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',1);
                        case 'MINI'
                            % [Aicelib,infoAiceall,valid_idx_ice] = crmsab_loadAicelibconv_mini(opticelib,basenameWA,c,bands_opt);
                        otherwise
                            error('Undefined OPT_LIBRARY %s',optLibrary);
                    end
                    ancillaries(c).Aicelib = Aicelib;
                    ancillaries(c).infoAicelib = infoAiceall(valid_idx_ice);
                else
                    switch upper(optLibrary)
                        case 'FULL'
                            [Asurficelib] = crmsab_loadAsurficelibconv(optsurficelib,basenameWA,c,...
                                bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',1);
                        case 'MINI'
                            % [Aicelib] = crmsab_loadAicelibconv_mini(opticelib,basenameWA,c,bands_opt);
                        otherwise
                            error('Undefined OPT_LIBRARY %s',optLibrary);
                    end
                end
            else
                Asurficelib = [];
            end
            if strcmpi(precision,'single')
                Alib = single(Alib); Aicelib = single(Aicelib);
            end
            if verLessThan('matlab','9.4')
                [ logt_est,logYifc_cor,logAB,logBg,logIce,logYifc_cor_ori,logYifc_isnan,Xt_c,Xlib_c,Xice_c,vldpxl_c]...
                    = sabcondc_v3l1_pub_R2017b(Alib,logYif(:,:,c),WAb(:,c),logT_extrap(:,:,c),GP(:,:,c),...
                          'LAMBDA_A',lambda_a,'NITER',nIter,'PRECISION',precision,'GPU',gpu,...
                          'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                          'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                          'Aicelib',Aicelib,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                          'THRESHOLD_BADSPC',th_badspc);

            else
            [ logt_est,logYifc_cor,logAB,logBg,logIce,logYifc_cor_ori,logYifc_isnan,Xt_c,Xlib_c,Xice_c,vldpxl_c]...
                = sabcondc_v3l1_pub(Alib,logYif(:,:,c),WAb(:,c),logT_extrap(:,:,c),GP(:,:,c),...
                  'LAMBDA_A',lambda_a,'NITER',nIter,'PRECISION',precision,'GPU',gpu,...
                  'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                  'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                  'Aicelib',Aicelib,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                  'THRESHOLD_BADSPC',th_badspc);
            end

            Yif_cor(lBool,c,bBool) = reshape(logYifc_cor',[nL,1,nB]);
            Yif_cor_ori(lBool,c,bBool) = reshape(logYifc_cor_ori',[nL,1,nB]);
            T_est(bBool,c) = logt_est;
            Bg_est(lBool,c,bBool) = reshape(logBg',[nL,1,nB]);
            AB_est(lBool,c,bBool) = reshape(logAB',[nL,1,nB]);
            Yif_isnan(lBool,c,bBool) = reshape(logYifc_isnan',[nL,1,nB]);
            if ~isempty(opticelib)
                Ice_est(lBool,c,bBool) = reshape(logIce',[nL,1,nB]);
            end
            ancillaries(c).Xt = Xt_c;
            ancillaries(c).Xlib = Xlib_c;
            ancillaries(c).Xice = Xice_c;
            Valid_pixels(lBool,c) = vldpxl_c';
            
            % toc;
        end

        Yif_cor = exp(Yif_cor); T_est = exp(T_est);
        Bg_est = exp(Bg_est); AB_est = exp(AB_est);
        Yif_cor_ori = exp(Yif_cor_ori); 
        if ~isempty(opticelib)
            Ice_est = exp(Ice_est);
        end
    
    case {'CPU_2','GPU_2','GPU_BATCH_2','CPU_3','GPU_3','GPU_BATCH_3','CPU_4','GPU_4','GPU_BATCH_4'}
        %%
        Yif_cor = nan([nBall,nLall,nCall],precision);
        Yif_isnan = nan([nBall,nLall,nCall]);
        switch upper(PROC_MODE)
            case {'CPU_2','GPU_2','GPU_BATCH_2'}
                T_est = nan([nBall,1,nCall],precision);
            case {'CPU_3','GPU_3','GPU_BATCH_3'}
                mc = nan([nBall,1,nCall],precision);
            case {'CPU_4','GPU_4','GPU_BATCH_4'}
                T_est = nan([nBall,1,nCall],precision);
                mc = nan([nBall,1,nCall],precision);
        end
        Bg_est = nan([nBall,nLall,nCall],precision);
        AB_est = nan([nBall,nLall,nCall],precision);
        badspcs = true([1,nLall,nCall]);
        bp_est_bools = false([nBall,1,nCall]);
        Ice_est = nan([nBall,nLall,nCall],precision);
        SurfIce_est = nan([nBall,nLall,nCall],precision);
        Yif_cor_ori = nan([nBall,nLall,nCall],precision);
        ancillaries = struct('Xt',cell(nCall,1),'Xlib',cell(nCall,1),...
                             'Xice',cell(nCall,1),...
                             'Alib',cell(nCall,1),'Aicelib',cell(nCall,1), 'Asurficelib',cell(nCall,1), ...
                             'infoAlib',cell(nCall,1),'infoAicelib',cell(nCall,1),'infoAsurficelib',cell(nCall,1), ...
                             'AlogT',cell(nCall,1));
        n_batch = ceil(length(Columns_valid)/batch_size);
        
        for ni = 1:n_batch
            tic;
            if ni~=n_batch
                Columns = Columns_valid((1+batch_size*(ni-1)):(batch_size*ni));
            elseif ni==n_batch
                Columns = Columns_valid((1+batch_size*(ni-1)):length(Columns_valid));
            end
            for i = 1:length(Columns)
                c = Columns(i);
                if Alib_out && (ni==1) && (i==1)
                    switch upper(optLibrary)
                        case 'FULL'
                            [Alib,infoAall,valid_idx] = crmsab_loadAlibconv(...
                                optLibs,basenameWA,optInterpid,c,...
                                bands_opt,WAb(:,c),cntRmvl);
                        case 'MINI'
                            [Alib,infoAall,valid_idx] = crmsab_loadAlibconv_mini( ...
                                optLibs,basenameWA,optInterpid,c,bands_opt);
                        otherwise
                            error('Undefined OPT_LIBRARY %s',optLibrary);
                    end
                    ancillaries(c).Alib     = Alib;
                    ancillaries(c).infoAlib = infoAall(valid_idx);
                else
                    switch upper(optLibrary)
                        case 'FULL'
                            [Alib] = crmsab_loadAlibconv(...
                                optLibs,basenameWA,optInterpid,c,...
                                bands_opt,WAb(:,c),cntRmvl);
                        case 'MINI'
                            [Alib] = crmsab_loadAlibconv_mini( ...
                                optLibs,basenameWA,optInterpid,c,bands_opt);
                        otherwise
                            error('Undefined OPT_LIBRARY %s',optLibrary);
                    end
                end
                if ~isempty(opticelib)
                    if Alib_out
                        switch upper(optLibrary)
                            case 'FULL'
                                [Aicelib,infoAiceall,valid_idx_ice]...
                                    = crmsab_loadAicelibconv(opticelib,basenameWA,c,...
                                    bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                            case 'MINI'
                                [Aicelib,infoAiceall,valid_idx_ice] ...
                                    = crmsab_loadAicelibconv_mini(opticelib,basenameWA,c,bands_opt);
                            otherwise
                                error('Undefined OPT_LIBRARY %s',optLibrary);
                        end
                        ancillaries(c).Aicelib     = Aicelib;
                        ancillaries(c).infoAicelib = infoAiceall(valid_idx_ice);
                    else
                        switch upper(optLibrary)
                            case 'FULL'
                                [Aicelib] = crmsab_loadAicelibconv(...
                                    opticelib,basenameWA,c,...
                                    bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',0);
                            case 'MINI'
                                [Aicelib] = crmsab_loadAicelibconv_mini(...
                                    opticelib,basenameWA,c,bands_opt);
                            otherwise
                                error('Undefined OPT_LIBRARY %s',optLibrary);
                        end
                    end
                else
                    Aicelib = [];
                end
                if ~isempty(optsurficelib)
                    if Alib_out
                        switch upper(optLibrary)
                            case 'FULL'
                                [Asurficelib,infoAsurficeall,valid_idx_surfice]...
                                    = crmsab_loadAsurficelibconv(optsurficelib,basenameWA,c,...
                                    bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',1);
                            case 'MINI'
                                % [Aicelib,infoAiceall,valid_idx_ice] ...
                                %     = crmsab_loadAicelibconv_mini(opticelib,basenameWA,c,bands_opt);
                            otherwise
                                error('Undefined OPT_LIBRARY %s',optLibrary);
                        end
                        ancillaries(c).Asurficelib     = Asurficelib;
                        ancillaries(c).infoAsurficelib = infoAsurficeall(valid_idx_surfice);
                    else
                        switch upper(optLibrary)
                            case 'FULL'
                                [Asurficelib] = crmsab_loadAsurficelibconv(...
                                    optsurficelib,basenameWA,c,...
                                    bands_opt,WAb(:,c),'overwrite',0,'CNTRMVL',1);
                            case 'MINI'
                                % [Aicelib] = crmsab_loadAicelibconv_mini(...
                                %     opticelib,basenameWA,c,bands_opt);
                            otherwise
                                error('Undefined OPT_LIBRARY %s',optLibrary);
                        end
                    end
                else
                    Asurficelib = [];
                end
                if i==1
                    Alibs = Alib;
                    Aicelibs = Aicelib;
                    Asurficelibs = Asurficelib;
                else
                    Alibs = cat(3,Alibs,Alib);
                    Aicelibs = cat(3,Aicelibs,Aicelib);
                    Asurficelibs = cat(3, Asurficelibs, Asurficelib);
                end
            end
            if strcmpi(precision,'single')
                Alibs = single(Alibs); Aicelibs = single(Aicelibs);
                Asurficelibs = single(Asurficelibs);
            end
            
            switch weight_mode
                case 0
                    if verLessThan('matlab','9.4')
                        [Yif_cor(bBool,lBool,Columns),T_est(bBool,1,Columns),...
                        AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns),...
                        Yif_isnan(bBool,lBool,Columns),Xt_c,Xlib_c,Xice_c,badspcs(1,lBool,Columns)]...
                            = sabcondc_v3l1_gpu_batch_tagv1_R2017b(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                                  logT_extrap(:,:,Columns),...
                                  BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                                  'Aicelib',Aicelibs,'nIter',nIter,...
                                  'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                                  'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                                  'gpu',gpu,'WEIGHT_MODE',weight_mode,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                                  'THRESHOLD_BADSPC',th_badspc);
                    else
                        if is_debug
                             [Yif_cor(bBool,lBool,Columns),T_est(bBool,1,Columns),...
                                AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns),...
                                Yif_isnan(bBool,lBool,Columns),Xt_c,Xlib_c,Xice_c,badspcs(1,lBool,Columns)]...
                            = sabcondc_v3l1_gpu_batch_tagv1_debug(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                                      logT_extrap(:,:,Columns),...
                                      BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                                      'Aicelib',Aicelibs,'nIter',nIter,...
                                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                                      'gpu',gpu,'WEIGHT_MODE',weight_mode,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                                      'THRESHOLD_BADSPC',th_badspc,'DEBUG',is_debug,'debug_l_plot',debug_l_plot);
                        else
                            [Yif_cor(bBool,lBool,Columns),T_est(bBool,1,Columns),...
                                AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns),...
                                Yif_isnan(bBool,lBool,Columns),Xt_c,Xlib_c,Xice_c,badspcs(1,lBool,Columns)]...
                            = sabcondc_v3l1_gpu_batch_tagv1(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                                      logT_extrap(:,:,Columns),...
                                      BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                                      'Aicelib',Aicelibs,'nIter',nIter,...
                                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                                      'gpu',gpu,'WEIGHT_MODE',weight_mode,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                                      'THRESHOLD_BADSPC',th_badspc);
                        end
                    end
                case 1
                    if verLessThan('matlab','9.4')
                        error('This version of MATLAB is not supported for Weight mode=2. Update at least to R2018a');
                    end
                    switch upper(PROC_MODE)
                         case {'CPU_4','GPU_4','GPU_BATCH_4'}
                             [Yif_cor(bBool,lBool,Columns),T_est(bBool,1,Columns),mc(bBool,1,Columns),...
                                AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns),...
                                Yif_isnan(bBool,lBool,Columns),Xt_c,Xlib_c,Xice_c,badspcs(1,lBool,Columns)]...
                            = sabcondc_v3l1_gpu_batch_mc2(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                                      logT_extrap(:,:,Columns),...
                                      BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                                      'Aicelib',Aicelibs,'nIter',nIter,...
                                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                                      'gpu',gpu,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                                      'THRESHOLD_BADSPC',th_badspc,...
                                      'WEIGHT_MODE',weight_mode,...
                                      'STDL1_IFDF',stdl1_ifdf(:,:,Columns),'SFIMG',SFimg(:,:,Columns),...
                                      'WA_UM_PITCH',WA_um_pitch(:,:,Columns),'LBL',TRRIFdata.lbl,'FFC_MODE',ffc_mode,...
                                      'DEBUG',is_debug,...
                                      'Bands_Bias_MAD',bands_bias_mad);
                        case {'CPU_3','GPU_3','GPU_BATCH_3'}
                             [Yif_cor(bBool,lBool,Columns),mc(bBool,1,Columns),...
                                AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns),...
                                Yif_isnan(bBool,lBool,Columns),Xt_c,Xlib_c,Xice_c,badspcs(1,lBool,Columns),bp_est_bools(bBool,1,Columns)]...
                            = sabcondc_v3l1_gpu_batch_mc(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                                      logT_extrap(:,:,Columns),...
                                      BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                                      'Aicelib',Aicelibs,'nIter',nIter,...
                                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                                      'gpu',gpu,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                                      'THRESHOLD_BADSPC',th_badspc,...
                                      'WEIGHT_MODE',weight_mode,...
                                      'STDL1_IFDF',stdl1_ifdf(:,:,Columns),'SFIMG',SFimg(:,:,Columns),...
                                      'WA_UM_PITCH',WA_um_pitch(:,:,Columns),'LBL',TRRIFdata.lbl,'FFC_MODE',ffc_mode,...
                                      'DEBUG',is_debug,...
                                      'Bands_Bias_MAD',bands_bias_mad);
                        case {'CPU_2','GPU_2','GPU_BATCH_2'}
                            [Yif_cor(bBool,lBool,Columns),T_est(bBool,1,Columns),...
                                AB_est(bBool,lBool,Columns),Bg_est(bBool,lBool,Columns),Ice_est(bBool,lBool,Columns), SurfIce_est(bBool,lBool,Columns),...
                                Yif_isnan(bBool,lBool,Columns),Xt_c,Xlib_c,Xice_c,Xsurfice_c,badspcs(1,lBool,Columns),bp_est_bools(bBool,1,Columns)]...
                            = sabcondc_v3l1_gpu_batch(logYif(:,:,Columns),WAb(:,Columns),Alibs,...
                                      logT_extrap(:,:,Columns),...
                                      BP(:,:,Columns),'lambda_a',lambda_a,'precision',precision,...
                                      'Aicelib',Aicelib, 'Asurficelib', Asurficelib,'nIter',nIter,...
                                      'lambda_a_ice', lambda_a_ice, 'lambda_a_surfice', lambda_a_surfice, ...
                                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb,...
                                      'gpu',gpu,'LAMBDA_UPDATE_RULE',lambda_update_rule,...
                                      'THRESHOLD_BADSPC',th_badspc,...
                                      'WEIGHT_MODE',weight_mode,...
                                      'STDL1_IFDF',stdl1_ifdf(:,:,Columns),'SFIMG',SFimg(:,:,Columns),...
                                      'WA_UM_PITCH',WA_um_pitch(:,:,Columns),'LBL',TRRIFdata.lbl,'FFC_MODE',ffc_mode,...
                                      'DEBUG',is_debug,...
                                      'Bands_Bias_MAD',bands_bias_mad,...
                                      'T_UPDATE',t_update,'LOGT_NEG',logT_neg,'logt_relax',logt_relax,...
                                      'opt_bands_ignore_init',opt_bands_ignore_init,'debug_l_plot',debug_l_plot, ...
                                      'Include_Ice', include_ice, 'THRESHOLD_SNR', thrd_snr);
                    end
            end
            switch upper(PROC_MODE)
                case {'GPU_BATCH_2'}
                    if include_ice
                        Yif_cor_ori(bBool,lBool,Columns)...
                            = logYif(:,:,Columns) ...
                            - gather(pagefun(@mtimes,gpuArray(T_est(bBool,1,Columns)),gpuArray(Xt_c)));
                    else
                        Yif_cor_ori(bBool,lBool,Columns)...
                            = logYif(:,:,Columns) ...
                            - gather(pagefun(@mtimes,gpuArray(T_est(bBool,1,Columns)),gpuArray(Xt_c)))...
                            -Ice_est(bBool,lBool,Columns);
                    end
                case {'GPU_BATCH_3'}
                    Yif_cor_ori(bBool,lBool,Columns)...
                    = logYif(:,:,Columns) ...
                    - gather(pagefun(@mtimes,gpuArray(logT_extrap(:,:,Columns)),gpuArray(Xt_c)))...
                    -Ice_est(bBool,lBool,Columns) - mc(bBool,:,Columns);
                case {'CPU_3','GPU_3'}
                    Yif_cor_ori(bBool,lBool,Columns) ...
                    = logYif(:,:,Columns)-logT_extrap(:,:,Columns)*Xt_c...
                      -Ice_est(bBool,lBool,Columns) - mc(bBool,:,Columns);
                case {'GPU_BATCH_4'}
                    Yif_cor_ori(bBool,lBool,Columns)...
                    = logYif(:,:,Columns) ...
                    - gather(pagefun(@mtimes,gpuArray(T_est(bBool,1,Columns)),gpuArray(Xt_c)))...
                    -Ice_est(bBool,lBool,Columns) - mc(bBool,:,Columns);
                case {'CPU_4','GPU_4'}
                    Yif_cor_ori(bBool,lBool,Columns) ...
                    = logYif(:,:,Columns)-T_est(bBool,1,Columns)*Xt_c...
                      -Ice_est(bBool,lBool,Columns) - mc(bBool,:,Columns);
                otherwise
                    if include_ice
                        Yif_cor_ori(bBool,lBool,Columns) ...
                        = logYif(:,:,Columns)-T_est(bBool,1,Columns)*Xt_c;
                    else
                        Yif_cor_ori(bBool,lBool,Columns) ...
                        = logYif(:,:,Columns)-T_est(bBool,1,Columns)*Xt_c...
                          -Ice_est(bBool,lBool,Columns);
                    end
            end
            
            Xt_c_cell   = squeeze(num2cell(Xt_c,[1,2]));
            Xlib_c_cell = squeeze(num2cell(Xlib_c,[1,2]));
            Xice_c_cell = squeeze(num2cell(Xice_c,[1,2]));
            Xsurfice_c_cell = squeeze(num2cell(Xsurfice_c,[1,2]));
            [ancillaries(Columns).Xt] = Xt_c_cell{:};
            [ancillaries(Columns).Xlib] = Xlib_c_cell{:};
            [ancillaries(Columns).Xice] = Xice_c_cell{:};
            [ancillaries(Columns).Xsurfice] = Xsurfice_c_cell{:};
            
            switch upper(PROC_MODE)
                case {'CPU_3','GPU_3','GPU_BATCH_3'}
                    logT_cell = squeeze(num2cell(logT_extrap(:,:,Columns),[1,2]));
                    [ancillaries(Columns).AlogT] = logT_cell{:};
                    
            end
            toc;
            
        end
        Yif_cor = permute(Yif_cor,[2,3,1]);
        Yif_cor_ori = permute(Yif_cor_ori,[2,3,1]);
        Yif_isnan = permute(Yif_isnan,[2,3,1]);
        Bg_est = permute(Bg_est,[2,3,1]);
        AB_est = permute(AB_est,[2,3,1]);
        Ice_est = permute(Ice_est,[2,3,1]);
        SurfIce_est = permute(SurfIce_est, [2,3,1]);
        switch upper(PROC_MODE)
            case {'CPU_2','GPU_2','GPU_BATCH_2'}
                T_est  = squeeze(T_est);
            case {'CPU_3','GPU_3','GPU_BATCH_3'}
                mc = squeeze(mc);
            case {'CPU_4','GPU_4','GPU_BATCH_4'}
                T_est  = squeeze(T_est);
                mc = squeeze(mc);
        end
        bp_est_bools = permute(bp_est_bools,[2,3,1]);
        badspcs = permute(badspcs,[2, 3, 1]);
        Valid_pixels = ~badspcs;
    
        Yif_cor = exp(Yif_cor); 
        Bg_est = exp(Bg_est); AB_est = exp(AB_est); 
        Ice_est = exp(Ice_est);
        SurfIce_est = exp(SurfIce_est);
        Yif_cor_ori = exp(Yif_cor_ori);
        switch upper(PROC_MODE)
            case {'CPU_2','GPU_2','GPU_BATCH_2'}
                T_est = exp(T_est);
            case {'CPU_3','GPU_3','GPU_BATCH_3'}
                mc = exp(mc);
            case {'CPU_4','GPU_4','GPU_BATCH_4'}
                T_est = exp(T_est);
                mc = exp(mc);
        end
    otherwise
        error('Undefined PROC_MODE=%s',PROC_MODE);
        
end


switch lower(optBP)
    case {'pri','all'}
        BP = (BP1nanfull==1);
        GP = (GP1nanfull==1);
    case 'none'
        BP = permute(BP,[2,3,1]);
        GP = permute(GP,[2,3,1]);
    otherwise 
        error('optBP=%s is not defined.',optBP);
end



tend = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
fprintf('finish procceing, current time is %s.\n',tend);
fprintf('Processing time is %s\n',tend-tstart);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a setting file.
fname_setting = [basename_cr '_settings.txt'];
if save_file
    fprintf('Now saving...\n');
    setting_file_path = joinPath(save_dir,fname_setting);
    fid = fopen(setting_file_path,'w');
else
    fid = 1; % standard output
end
hostname = char(java.net.InetAddress.getLocalHost.getHostName);
dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
fprintf(fid,'Hostname: %s\n',hostname);
fprintf(fid,'Username: %s\n',username);
fprintf(fid,'Time finished: %s\n',dt);

% ## I/O OPTIONS #---------------------------------------------------------
fprintf(fid,'SAVE_FILE: %d\n',save_file);
fprintf(fid,'STORAGE_SAVING_LEVEL %s\n',storage_saving_level);
fprintf(fid,'SAVE_PDIR: %s\n',save_pdir);
fprintf(fid,'SAVE_DIR_YYYY_DOY: %d\n',save_dir_yyyy_doy);
fprintf(fid,'FORCE: %d\n',force);
fprintf(fid,'SKIP_IFEXIST: %d\n',skip_ifexist);
fprintf(fid, 'SUFFIX: %s\n', suffix);
fprintf(fid,'ADDITIONAL_SUFFIX: %s\n',additional_suffix);
fprintf(fid,'INTERLEAVE_OUT: %s\n',interleave_out);
fprintf(fid,'SUBSET_COLUMNS_OUT: %d\n', subset_columns_out);
fprintf(fid,'ALIB_OUT: %d\n',Alib_out);
fprintf(fid,'CROP_BANDS: %d\n',do_crop_bands);

% ## INPUT IMAGE OPTIONS #-------------------------------------------------
fprintf(fid,'OPT_IMG: %s\n',opt_img);
fprintf(fid,'IMG_CUBE is empty: %s',img_cube_isempty);
fprintf(fid,'IMG_CUBE_BAND_INVERSE: %d',img_cube_band_inverse);
fprintf(fid,'TRRY_PDIR: %s\n',dir_TRRX);
fprintf(fid,'FFC_IF_COUNTER: %d\n',ffc_counter);
fprintf(fid,'OBS_COUNTER: %d\n', obs_counter);
fprintf(fid,'OBS_COUNTER_DF: %d\n', obs_counter_df);

% ## GENERAL SABCOND OPTIONS #---------------------------------------------
fprintf(fid,'BANDS_OPT: %d\n',bands_opt);
fprintf(fid,'LINES:'); fprintf(fid,' %d',line_idxes); fprintf(fid, '\n');
fprintf(fid,'COLUMNS:'); fprintf(fid,' %d',column_idxes); fprintf(fid, '\n');
fprintf(fid,'METHODTYPE: %s\n',mt);
fprintf(fid,'OPTBP: %s\n',optBP);
fprintf(fid,'VERBOSE: %d\n', verbose);
fprintf(fid,'COLUMN_SKIP: %d\n',column_skip);
fprintf(fid,'WEIGHT_MODE: %d\n',weight_mode);
fprintf(fid,'LAMBDA_UPDATE_RULE: %s\n',lambda_update_rule);
fprintf(fid,'THRESHOLD_BADSPC: %f\n',th_badspc);
fprintf(fid,'FFC_MODE: %d\n',ffc_mode);
fprintf(fid,'OPT_BANDS_IGNORE_INIT: %s\n',opt_bands_ignore_init);
fprintf(fid,'INCLUDE_ICE: %d\n', include_ice);

% ## PRE-PROCESSING OPTIONS #----------------------------------------------
fprintf(fid,'CAL_BIAS_COR: %d\n', cal_bias_cor);

% ## TRANSMISSION SPECTRUM OPTIONS #---------------------------------------
fprintf(fid,'T_MODE: %d\n', t_mode);
fprintf(fid,'OBS_ID_T: %s\n',obs_id_T);
fprintf(fid,'VARARGIN_T:');
for i=1:length(varargin_T)
    if isnumeric(varargin_T{i})
        string_varargin_T = num2str(varargin_T{i});
    else
        string_varargin_T = varargin_T{i};
    end
    fprintf(fid,' "%s"', string_varargin_T);
    if i<length(varargin_T)
        fprintf(',');
    end
end
fprintf(fid, '\n');

% ## LIBRARY OPTIONS #-----------------------------------------------------
fprintf(fid,'OPT_LIBRARY: %d\n',optLibrary);
fprintf(fid,'CNTRMVL: %d\n',cntRmvl);
fprintf(fid,'OPTINTERPID: %d\n',optInterpid);
fprintf(fid,'OPT_CRISMSPCLIB: %d\n',optCRISMspclib);
fprintf(fid,'OPT_RELAB: %d\n',optRELAB);
fprintf(fid,'OPT_SPLIBUSGS: %d\n',optUSGSsplib);
fprintf(fid,'OPT_CRISMTYPELIB: %d\n',optCRISMTypeLib);
fprintf(fid,'OPT_ICELIB: %d\n',opticelib);
fprintf(fid,'OPT_SURFICELIB: %d\n',optsurficelib);

% ## SABCONDC OPTIONS #----------------------------------------------------
fprintf(fid,'NITER: %d\n',nIter);
fprintf(fid,'LAMBDA_A:'); fprintf(fid,' %f',lambda_a); fprintf(fid,'\n');
fprintf(fid,'LAMBDA_A_ICE:'); fprintf(fid,' %f',lambda_a_ice); fprintf(fid,'\n');
fprintf(fid,'LAMBDA_A_SURFICE:'); fprintf(fid,' %f',lambda_a_surfice); fprintf(fid,'\n');
fprintf(fid,'T_UPDATE: %d\n',t_update);
fprintf(fid,'LOGT_NEG: %d\n',logT_neg);
fprintf(fid,'LOGT_RELAX: %d\n',logt_relax);
fprintf(fid,'THRD_SNR: %f\n',thrd_snr);

% ## PROCESSING OPTIONS #--------------------------------------------------
fprintf(fid,'PRECISION: %s\n',precision);
fprintf(fid,'PROC_MODE: %s\n',PROC_MODE);
fprintf(fid,'DEBUG: %d\n',is_debug);

% ## ETCETERA #------------------------------------------------------------
% fprintf(fid,'GAUSSSIGMA: %f\n',gausssigma);

if fid>1
    fclose(fid);
    chmod777(setting_file_path,verbose);
end

%% dump to settings
settings = [];
settings.hostname = hostname;
settings.username = username;
settings.dt = dt;
% ## I/O OPTIONS #---------------------------------------------------------
settings.save_file = save_file;
settings.storage_saving_level = storage_saving_level;
settings.save_pdir = save_pdir;
settings.save_dir_yyyy_doy = save_dir_yyyy_doy;
settings.force = force;
settings.skip_ifexist = skip_ifexist;
settings.suffix = suffix;
settings.additional_suffix = additional_suffix;
settings.interleave_out = interleave_out;
settings.subset_columns_out = subset_columns_out;
settings.Alib_out = Alib_out;
settings.crop_bands = do_crop_bands;
% ## INPUT IMAGE OPTIONS #-------------------------------------------------
settings.opt_img = opt_img;
settings.img_cube_isempty = img_cube_isempty;
settings.img_cube_band_inverse = img_cube_band_inverse;
settings.trry_pdir = dir_TRRX;
settings.ffc_if_counter = ffc_counter;
settings.obs_counter = obs_counter;
settings.obs_counter_df = obs_counter_df;
% ## GENERAL SABCOND OPTIONS #---------------------------------------------
settings.bands_opt = bands_opt;
settings.lines = line_idxes;
settings.columns = column_idxes;
settings.mt = mt;
settings.optBP = optBP;
settings.verbose = verbose;
settings.column_skip = column_skip;
settings.weight_mode = weight_mode;
settings.lambda_update_rule = lambda_update_rule;
settings.th_badspc = th_badspc;
settings.ffc_mode = ffc_mode;
settings.opt_bands_ignore_init = opt_bands_ignore_init;
settings.include_ice = include_ice;
% ## PRE-PROCESSING OPTIONS #----------------------------------------------
settings.cal_bias_cor = cal_bias_cor;
% ## TRANSMISSION SPECTRUM OPTIONS #---------------------------------------
settings.t_mode = t_mode;
settings.obs_id_T = obs_id_T;
settings.varargin_T = varargin_T;
% ## LIBRARY OPTIONS #-----------------------------------------------------
settings.optLibrary = optLibrary;
settings.cntRmvl = cntRmvl;
settings.optInterpid = optInterpid;
settings.optCRISMspclib = optCRISMspclib;
settings.optRELAB = optRELAB;
settings.optUSGSsplib = optUSGSsplib;
settings.optCRISMTypeLib = optCRISMTypeLib;
settings.opticelib = opticelib;
settings.optsurficelib = optsurficelib;
% ## SABCONDC OPTIONS #----------------------------------------------------
settings.nIter = nIter;
settings.lambda_a = lambda_a;
settings.lambda_a_ice = lambda_a_ice;
settings.lambda_a_surfice = lambda_a_surfice;
settings.t_update = t_update;
settings.logT_neg = logT_neg;
settings.logt_relax = logt_relax;
settings.thrd_snr = thrd_snr;
% ## PROCESSING OPTIONS #--------------------------------------------------
settings.precision = precision;
settings.PROC_MODE = PROC_MODE;
settings.debug = is_debug;

%% Post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performing interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yif_cor_nr = Yif_cor;
if include_ice
    Yif_mdl = AB_est .* Bg_est .* SurfIce_est .* Ice_est;
else
    Yif_mdl = AB_est .* Bg_est .* SurfIce_est;
end
nan_cells = isnan(Yif_cor);
Yif_cor_nr(nan_cells) = Yif_mdl(nan_cells);

Valid_pixels_good = double(Valid_pixels);
Valid_pixels_good(Valid_pixels==0) = nan;
for bi=1:nBall
    Yif_cor_nr(:,:,bi) = Yif_cor_nr(:,:,bi) .* Valid_pixels_good;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create header files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hdr, mimics CAT file production
% first construct ENVI header files
hdr_cr = crism_const_cathdr(TRRIFdata,true,'DATE_TIME',dt);
hdr_cr.cat_history = suffix;
% update bbl
bbl = false(1,hdr_cr.bands);
bbl(bands) = true;
hdr_cr.bbl = bbl;

switch opt_img
    case 'if'
        hdr_cr.cat_input_files = [TRRIFdata.basename '.IMG'];
    case 'ra_if'
        hdr_cr.cat_input_files = [TRRRAIFdata.basename '_IF.IMG'];
    case {'TRRY','TRRC','TRRB','TRRD','TRRE'}
        hdr_cr.cat_input_files = [basenameTRRY '.IMG'];
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desmiling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desmiling is performed by simple linear interpolation.
% Only performed for continuous data over the bands.
fprintf('Performing de-smiling ...\n');
wa = WAdata.img;
[Yif_nr_ds] = crism_smile_correction(Yif_cor_nr,wa,hdr_cr.wavelength(:)*1000,bands);
[Yifmdl_ds] = crism_smile_correction(Yif_mdl   ,wa,hdr_cr.wavelength(:)*1000,bands);
[AB_est_ds] = crism_smile_correction(AB_est    ,wa,hdr_cr.wavelength(:)*1000,bands);
[Bg_est_ds] = crism_smile_correction(Bg_est    ,wa,hdr_cr.wavelength(:)*1000,bands);
[Ice_est_ds] = crism_smile_correction(Ice_est  ,wa,hdr_cr.wavelength(:)*1000,bands);
[SurfIce_est_ds] = crism_smile_correction(SurfIce_est  ,wa, hdr_cr.wavelength(:)*1000,bands);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crop bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if do_crop_bands
    Yif_cor     = Yif_cor(:,:,bands);
    Yif_cor_ori = Yif_cor_ori(:,:,bands);
    AB_est      = AB_est(:,:,bands);
    Bg_est      = Bg_est(:,:,bands);
    if ~isempty(opticelib)
        Ice_est = Ice_est(:,:,bands);
    end
    if ~isempty(optsurficelib)
        SurfIce_est = SurfIce_est(:,:,bands);
    end
    Yif_nr_ds   = Yif_nr_ds(:,:,bands);
    Yifmdl_ds   = Yifmdl_ds(:,:,bands);
    AB_est_ds   = AB_est_ds(:,:,bands);
    Bg_est_ds   = Bg_est_ds(:,:,bands);
    hdr_cr.wavelength = hdr_cr.wavelength(bands);
    hdr_cr.fwhm  = hdr_cr.fwhm(bands);
    hdr_cr.bbl   = hdr_cr.bbl(bands);
    hdr_cr.bands = length(bands);
    wa = wa(:,:,bands);
    BP = BP(:,:,bands);
    GP = GP(:,:,bands);
    hdr_cr.band_names = arrayfun(@(x) sprintf('Band %3d of %s',x,TRRIFdata.basename),(nBall-bands+1),...
    'UniformOutput',false);
end

%% SAVING OPERATIONS
if save_file
    switch upper(storage_saving_level)
        case 'NORMAL'
            save_output_img(single(Yif_cor),hdr_cr,basename_cr,save_dir,verbose);
            % fpath_hdr = joinPath(save_dir, [basename_cr '.hdr']);
            % fprintf('Saving %s ...\n',fpath_hdr);
            % envihdrwritex(hdr_cr,fpath_hdr,'OPT_CMOUT',false);
            % chmod777(fpath_hdr,verbose);
            % fprintf('Done\n');
            % fpath_img = joinPath(save_dir, [basename_cr '.img']);
            % fprintf('Saving %s ...\n',fpath_img);
            % envidatawrite(single(Yif_cor),fpath_img,hdr_cr);
            % chmod777(fpath_img,verbose);
            % fprintf('Done\n');
            
            basename_ori = [basename_cr '_ori'];
            save_output_img(single(Yif_cor_ori),hdr_cr,basename_ori,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_ori '.hdr']));
            % envihdrwritex(hdr_cr,joinPath(save_dir,[basename_ori '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_ori '.img']));
            % envidatawrite(single(Yif_cor_ori),joinPath(save_dir,[basename_ori '.img']),hdr_cr);
            % fprintf('Done\n');

            fpath_supple = joinPath(save_dir,[basename_cr '.mat']);
            wa = squeeze(wa)';
            fprintf('Saving %s ...\n',fpath_supple);

            switch upper(PROC_MODE)
                case {'CPU_3','GPU_3','GPU_BATCH_3'}
                    save(fpath_supple,'wa','bands','line_idxes','mc','BP','GP','bp_est_bools',...
                        'ancillaries','Valid_pixels','settings');
                case {'CPU_4','GPU_4','GPU_BATCH_4'}
                    save(fpath_supple,'wa','bands','line_idxes','T_est','mc','BP','GP','bp_est_bools',...
                        'ancillaries','Valid_pixels','settings');
                otherwise
                    save(fpath_supple,'wa','bands','line_idxes','T_est','BP','GP','bp_est_bools',...
                        'ancillaries','Valid_pixels','settings');
            end
            chmod777(fpath_supple,verbose);
            fprintf('Done\n');

            basename_Bg = [basename_cr '_Bg'];
            save_output_img(single(Bg_est),hdr_cr,basename_Bg,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg '.hdr']));
            % envihdrwritex(hdr_cr,joinPath(save_dir,[basename_Bg '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg '.img']));
            % envidatawrite(single(Bg_est),joinPath(save_dir, [basename_Bg '.img']),hdr_cr);
            % fprintf('Done\n');

            basename_AB = [basename_cr '_AB'];
            save_output_img(single(AB_est),hdr_cr,basename_AB,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB '.hdr']));
            % envihdrwritex(hdr_cr,joinPath(save_dir, [basename_AB '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB '.img']));
            % envidatawrite(single(AB_est),joinPath(save_dir, [basename_AB '.img']),hdr_cr);
            % fprintf('Done\n');

            if ~isempty(opticelib)
                basename_Ice = [basename_cr '_Ice'];
                save_output_img(single(Ice_est),hdr_cr,basename_Ice,save_dir,verbose);
                % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Ice '.hdr']));
                % envihdrwritex(hdr_cr,joinPath(save_dir, [basename_Ice '.hdr']),'OPT_CMOUT',false);
                % fprintf('Done\n');
                % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Ice '.img']));
                % envidatawrite(single(Ice_est),joinPath(save_dir, [basename_Ice '.img']),hdr_cr);
                % fprintf('Done\n');
                basename_Ice_ds = [basename_Ice  '_ds'];
                hdr_ds = hdr_cr;
                dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
                hdr_ds.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
                hdr_ds.cat_history = [hdr_cr.cat_history '_ds'];
                save_output_img(single(Ice_est_ds), hdr_ds, ...
                    basename_Ice_ds, save_dir, verbose);
            end

            if ~isempty(optsurficelib)
                basename_SurfIce = [basename_cr '_SurfIce'];
                save_output_img(single(SurfIce_est),hdr_cr,basename_SurfIce,save_dir,verbose);
                % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Ice '.hdr']));
                % envihdrwritex(hdr_cr,joinPath(save_dir, [basename_Ice '.hdr']),'OPT_CMOUT',false);
                % fprintf('Done\n');
                % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Ice '.img']));
                % envidatawrite(single(Ice_est),joinPath(save_dir, [basename_Ice '.img']),hdr_cr);
                % fprintf('Done\n');
                basename_SurfIce_ds = [basename_SurfIce  '_ds'];
                hdr_ds = hdr_cr;
                dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
                hdr_ds.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
                hdr_ds.cat_history = [hdr_cr.cat_history '_ds'];
                save_output_img(single(SurfIce_est_ds), hdr_ds, ...
                    basename_SurfIce_ds, save_dir,verbose);
            end
            
            % Nan Replaced data (nr)
            basename_cr_nr = [basename_cr '_nr'];
            hdr_cr_nr = hdr_cr;
            dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
            hdr_cr_nr.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
            hdr_cr_nr.cat_history = [hdr_cr_nr.cat_history '_nr'];
            save_output_img(single(Yif_cor_nr),hdr_cr_nr,basename_cr_nr,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr '.hdr']));
            % envihdrwritex(hdr_cr_nr,joinPath(save_dir,[basename_cr_nr '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr '.img']));
            % envidatawrite(single(Yif_cor_nr),joinPath(save_dir,[basename_cr_nr '.img']),hdr_cr);
            % fprintf('Done\n');

            basename_nr_ds = [basename_cr_nr '_ds'];
            hdr_nr_ds = hdr_cr_nr;
            dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
            hdr_nr_ds.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
            hdr_nr_ds.cat_history = [hdr_cr_nr.cat_history '_ds'];
            save_output_img(single(Yif_nr_ds),hdr_nr_ds,basename_nr_ds,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_nr_ds '.hdr']));
            % envihdrwritex(hdr_nr_ds,joinPath(save_dir,[basename_nr_ds '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_nr_ds '.img']));
            % envidatawrite(single(Yif_nr_ds),joinPath(save_dir,[basename_nr_ds '.img']),hdr_cr);
            % fprintf('Done\n');
            
            basename_Bg_ds = [basename_Bg '_ds'];
            save_output_img(single(Bg_est_ds),hdr_nr_ds,basename_Bg_ds,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg_ds '.hdr']));
            % envihdrwritex(hdr_nr_ds,joinPath(save_dir,[basename_Bg_ds '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg_ds '.img']));
            % envidatawrite(single(Bg_est_ds),joinPath(save_dir, [basename_Bg_ds '.img']),hdr_cr);
            % fprintf('Done\n');

            basename_AB_ds = [basename_AB '_ds'];
            save_output_img(single(AB_est_ds),hdr_nr_ds,basename_AB_ds,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB_ds '.hdr']));
            % envihdrwritex(hdr_nr_ds,joinPath(save_dir, [basename_AB_ds '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB_ds '.img']));
            % envidatawrite(single(AB_est_ds),joinPath(save_dir, [basename_AB_ds '.img']),hdr_cr);
            % fprintf('Done\n');
            
        case 'HIGHEST'
            % With STORAGE_SAVING_LEVEL=HIGHEST, only desmiled nan replaced
            % image and models are saved.
            basename_nr_ds = [basename_cr '_nr_ds'];
            hdr_nr_ds = hdr_cr;
            dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
            hdr_nr_ds.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
            hdr_nr_ds.cat_history = [hdr_cr.cat_history '_nr_ds'];
            save_output_img(single(Yif_nr_ds),hdr_nr_ds,basename_nr_ds,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_nr_ds '.hdr']));
            % envihdrwritex(hdr_nr_ds,joinPath(save_dir,[basename_nr_ds '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_nr_ds '.img']));
            % envidatawrite(single(Yif_nr_ds),joinPath(save_dir,[basename_nr_ds '.img']),hdr_cr);
            % fprintf('Done\n');
            
            basename_mdl_ds = [basename_cr '_mdl_ds'];
            hdr_mdl_ds = hdr_cr;
            dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
            hdr_mdl_ds.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
            hdr_mdl_ds.cat_history = [hdr_cr.cat_history '_mdl_ds'];
            save_output_img(single(Yifmdl_ds),hdr_mdl_ds,basename_mdl_ds,save_dir,verbose);
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_mdl_ds '.hdr']));
            % envihdrwritex(hdr_nr_ds,joinPath(save_dir,[basename_mdl_ds '.hdr']),'OPT_CMOUT',false);
            % fprintf('Done\n');
            % fprintf('Saving %s ...\n',joinPath(save_dir, [basename_mdl_ds '.img']));
            % envidatawrite(single(Yifmdl_ds),joinPath(save_dir, [basename_mdl_ds '.img']),hdr_cr);
            % fprintf('Done\n');
    end
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

if save_file
    diary off;
    chmod777(fpath_diary,verbose);
    % if isunix
    %     system(['chmod 777 ' fpath_diary]);
    %     if verbose, fprintf('"%s": permission is set to 777.\n',fpath_diary); end
    % end
end

%% Construct output
out = [];
if nargout==0
    
elseif nargout==1
    out.obs_id = obs_id;
    out.TRR3IF = TRRIFdata.basename;
    
    WA = squeeze(wa)';
    % take the subset of the columns
    if subset_columns_out
        Yif_cor     = Yif_cor(:,column_idxes,:);
        Yif_cor_nr  = Yif_cor_nr(:,column_idxes,:);
        Yif_cor_ori = Yif_cor_ori(:,column_idxes,:);
        Yif_nr_ds   = Yif_nr_ds(:,column_idxes,:);
        Yif_isnan   = Yif_isnan(:,column_idxes,:);
        AB_est      = AB_est(:,column_idxes,:);
        Bg_est      = Bg_est(:,column_idxes,:);
        AB_est_ds   = AB_est_ds(:,column_idxes,:);
        Bg_est_ds   = Bg_est_ds(:,column_idxes,:);
        Yifmdl_ds   = Yifmdl_ds(:,column_idxes,:);
        WA          = WA(:,column_idxes);
        switch upper(PROC_MODE)
            case {'CPU_3','GPU_3','GPU_BATCH_3'}
                mc = mc(:,column_idxes);
            case {'CPU_4','GPU_4','GPU_BATCH_4'}
                mc = mc(:,column_idxes);
                T_est = T_est(:,column_idxes);
            otherwise
                T_est = T_est(:,column_idxes);
        end
        Valid_pixels= Valid_pixels(:,column_idxes);
        GP = GP(:,column_idxes,:);
        BP = BP(:,column_idxes,:);
        bp_est_bools = bp_est_bools(:,column_idxes,:);
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
    Yif_nr_ds   = permute(Yif_nr_ds  ,prmt_ordr);
    Yif_isnan   = permute(Yif_isnan,  prmt_ordr);
    AB_est      = permute(AB_est,     prmt_ordr);
    Bg_est      = permute(Bg_est,     prmt_ordr);
    AB_est_ds   = permute(AB_est_ds,  prmt_ordr);
    Bg_est_ds   = permute(Bg_est_ds,  prmt_ordr);
    Yifmdl_ds   = permute(Yifmdl_ds,  prmt_ordr);
    Valid_pixels= permute(Valid_pixels,prmt_ordr);
    GP = permute(GP,prmt_ordr);
    BP = permute(BP,prmt_ordr);
    bp_est_bools = permute(bp_est_bools,prmt_ordr);
    if ~isempty(opticelib)
        Ice_est = permute(Ice_est,    prmt_ordr);
    end
    if ~isempty(optsurficelib)
        SurfIce_est = permute(SurfIce_est,    prmt_ordr);
    end

    out.Yif_cor     = Yif_cor;
    out.Yif_cor_nr  = Yif_cor_nr;
    out.Yif_cor_ori = Yif_cor_ori;
    out.Yif_nr_ds   = Yif_nr_ds;
    out.AB_est      = AB_est;
    out.Bg_est      = Bg_est;
    out.AB_est_ds   = AB_est_ds;
    out.Bg_est_ds   = Bg_est_ds;
    out.Yifmdl_ds   = Yifmdl_ds;
    if ~isempty(opticelib)
        out.Ice_est = Ice_est;
    end
    if ~isempty(opticelib)
        out.SurfIce_est = SurfIce_est;
    end
    switch upper(PROC_MODE)
        case {'CPU_3','GPU_3','GPU_BATCH_3'}
            out.mc = mc;
        case {'CPU_4','GPU_4','GPU_BATCH_4'}
            out.mc = mc;
            out.T_est = T_est;
        otherwise
            out.T_est = T_est;
    end
    out.bp_est_bools = bp_est_bools;
    out.Yif_isnan    = Yif_isnan;
    out.Valid_pixels = Valid_pixels;
    out.WA           = WA;
    out.WV           = hdr_cr.wavelength;
    out.lines        = line_idxes;
    out.columns      = column_idxes;
    out.bands        = bands;
    out.interleave_out     = interleave_out;
    out.subset_columns_out = subset_columns_out;
    out.GP = GP;
    out.BP = BP;
    
    if subset_columns_out
        ancillaries  = ancillaries(column_idxes);
    end
    out.ancillaries  = ancillaries;
    out.settings = settings;
    
end


end

function save_output_img(img,hdr,basename,save_dir,verbose)
fpath_hdr = joinPath(save_dir, [basename '.hdr']);
fprintf('Saving %s ...\n',fpath_hdr);
envihdrwritex(hdr,fpath_hdr,'OPT_CMOUT',false);
chmod777(fpath_hdr,verbose);
fprintf('Done\n');
fpath_img = joinPath(save_dir, [basename '.img']);
fprintf('Saving %s ...\n',fpath_img);
envidatawrite(img,fpath_img,hdr);
chmod777(fpath_img,verbose);
fprintf('Done\n');
end

