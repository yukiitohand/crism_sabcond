function [Yif_cor,Bg_est,AB_est,T_est] = sabcondv3_pub(obs_id,varargin)
% [Yif_cor,Bg_est,AB_est,T_est] = sabcondv3_pub(obs_id,varargin)
%   This function performs simultaneous de-noising and atmospheric
%   correction of CRISM data. Currently 

% INPUT parameters

% huwacb option
nIter = 5;
lambda_a = 0.01;
% band option
opt_img = 'TRRB';
bands_opt = 4;
% library options
optCRISMspclib = 1;
optRELAB = 1;
optUSGSsplib = 1;
optCRISMTypeLib = 2;
opticelib = '';
cntRmvl = 1;
optInterpid = 1;
t_mode = 2;
% file name option
mt = 'sabcondv3'; % method type
additional_suffix = 'v1';
save_pdir = './resu/';
save_dir_yyyy_doy = false;
force = false;
gausssigma = 0.6;
optBP = 'pri'; %{'pri','all','none'}
skip_ifexist = false;

% ffc counter
ffc_counter = 1;
% specifying lines to be processed
lls = [];

precision = 'double';

gpu = false;
verbose = 0;

global crism_env_vars
dir_yuk = crism_env_vars.dir_YUK;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'FFC_IF_COUNTER'
                ffc_counter = varargin{i+1};
            case 'NITER'
                nIter = varargin{i+1};
            case 'LAMBDA_A'
                lambda_a = varargin{i+1};
            case 'BANDS_OPT'
                bands_opt = varargin{i+1};
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
            case 'OPT_IMG'
                opt_img = varargin{i+1};
            case 'CNTRMVL'
                cntRmvl = varargin{i+1};
            case 'T_MODE'
                t_mode = varargin{i+1};
            case 'METHODTYPE'
                mt = varargin{i+1};
            case 'ADDITIONAL_SUFFIX'
                additional_suffix = varargin{i+1};
            case 'OPTINTERPID'
                optInterpid = varargin{i+1};
            case 'SAVE_PDIR'
                save_pdir = varargin{i+1};
            case 'SAVE_DIR_YYYY_DOY'
                save_dir_yyyy_doy = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case 'OPTBP'
                optBP = varargin{i+1};
            case 'SKIP_IFEXIST'
                skip_ifexist = varargin{i+1};
            case 'TRRY_PDIR'
                dir_yuk = varargin{i+1};
            case 'LLS'
                lls = varargin{i+1};
            case 'PRECISION'
                precision = varargin{i+1};
            case 'GPU'
                gpu = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

if force && skip_ifexist
    error('You are forcing or skipping? Not sure what you want');
end

if ~exist(save_pdir,'dir')
    mkdir(save_pdir);
end

fprintf('lambda_a:%f\n',lambda_a);
fprintf('nIter:%d\n',nIter);
fprintf('t_mode:%d\n',t_mode);
fprintf('gauss_sigma:%f\n',gausssigma);
fprintf('opt_img: %s\n',opt_img);
fprintf('optBP: %s\n',optBP);
fprintf('opt_icelib: %d\n',opticelib);

bands = genBands(bands_opt);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,'');

% open log file
username = char(java.lang.System.getProperty('user.name'));
fname = sprintf('log_%s_%s.txt',username,datetime('now','TimeZone','local','Format','yyyyMMdd'));
diary(joinPath(save_pdir,fname));

fprintf('Current directory:%s\n',pwd);

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

[suffix] = const_suffix_v2(mt,cntRmvl,libprefix,optInterpid,bands_opt,nIter,additional_suffix);
fprintf('suffix will be \n%s.\n',suffix);

switch opt_img
    case 'if'
        basename_cr = [crism_obs.info.basenameIF suffix];
    case 'ra_if'
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


if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% Read image and ancillary data and format them for processing
TRRIFdata.load_basenamesCDR();
WAdata = TRRIFdata.readCDR('WA'); WAdata.readimgi();
nLall = TRRIFdata.hdr.lines; nCall = TRRIFdata.hdr.samples; nBall = TRRIFdata.hdr.bands;
if isempty(lls), lls = 1:nLall; end
nL = length(lls); nB = length(bands);
lBool = false(nLall,1); lBool(lls) = true;
bBool = false(nBall,1); bBool(bands) = true;

basenameWA = WAdata.basename;
WA = squeeze(WAdata.img(:,:,bands))';

switch opt_img
    case 'if'
        if TRRIF_is_empty
            error('TRR I/F does not exist.');
        else
            Yif = TRRIFdata.readimgi();
        end
    case 'ra_if'
        TRRRAIFdata = CRISMdataCAT([crism_obs.info.basenameRA '_IF'],crism_obs.info.dir_trdr,'ra_if');
        Yif = TRRRAIFdata.readimgi();
    case {'TRRY','TRRB','TRRC'}
        d_IoF = joinPath(dir_yuk, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        TRRYIFdata = CRISMdata(basenameTRRY,d_IoF);
        Yif = TRRYIFdata.readimgi();
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

Yif = Yif(lls,:,bands);
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
    WA     = single(WA);
end

% evaluate valid columns
valid_columns = ~all(isnan(WA),1);
Columns = 1:nCall;
Columns_valid = Columns(valid_columns);

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

switch PROC_MODE
    case 0
        Yif_cor = nan([nLall,nCall,nBall],precision);
        T_est = nan([nBall,nCall],precision);
        Bg_est = nan([nLall,nCall,nBall],precision);
        AB_est = nan([nLall,nCall,nBall],precision);
        ancillaries = struct('X',cell(1,nCall),'lambda',cell(1,nCall),...
            'nIter',cell(1,nCall),'huwacb_func',cell(1,nCall),...
            'maxiter_huwacb',cell(1,nCall),'tol_huwacb',cell(1,nCall),...
            'gp_bool',cell(1,nCall));
        Yif_cor_ori = nan([nLall,nCall,nBall],precision);
        Valid_pixels = false([nLall,nCall]);
        for c = Columns_valid
            tic;
            [Alib] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,...
                bands_opt,WA(:,c),cntRmvl);
            if ~isempty(opticelib)
                [Aicelib] = loadlibc_crism_icelib(opticelib,basenameWA,c,...
                    bands_opt,WA(:,c),'overwrite',0,'CNTRMVL',0);
            else
                Aicelib = [];
            end
            if strcmpi(precision,'single')
                Alib = single(Alib); Aicelib = single(Aicelib);
            end

            [ logt_est,logYifc_cor,logAB,logBg,logIce,logYifc_cor_ori,~,ancillary,~,vldpxl_c]...
                = sabcondc_v3l1_pub(Alib,logYif(:,:,c),WA(:,c),logT_extrap(:,:,c),...
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
            ancillaries(c) = ancillary;
            Valid_pixels(lBool,c) = vldpxl_c';

            toc;
        end

        Yif_cor = exp(Yif_cor); T_est = exp(T_est);
        Bg_est = exp(Bg_est); AB_est = exp(AB_est);
        Yif_cor_ori = exp(Yif_cor_ori); logIce = exp(logIce);
        
    case 1
        Yif_cor = nan([nBall,nLall,nCall],precision);
        Yif_isnan = nan([nBall,nLall,nCall]);
        T_est = nan([nBall,1,nCall],precision);
        Bg_est = nan([nBall,nLall,nCall],precision);
        AB_est = nan([nBall,nLall,nCall],precision);
        Yif_cor_ori = nan([nBall,nLall,nCall],precision);
        X = [];

        for ni = 1:n_batch
            if ni~=n_batch
                cList = Columns_valid((1+batch_size*(ni-1)):(batch_size*ni));
            elseif ni==n_batch
                cList = Columns_valid((1+batch_size*(ni-1)):length(Columns_valid));
            end
            for i = 1:length(cList)
                c = cList(i);
                [Alib] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,...
                    bands_opt,WA(:,c),cntRmvl);
                if ~isempty(opticelib)
                [Aicelib] = loadlibc_crism_icelib(opticelib,basenameWA,c,...
                    bands_opt,WA(:,c),'overwrite',0,'CNTRMVL',0);
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
            [Yif_cor(bBool,lBool,cList),T_est(bBool,1,cList),...
                AB_est(bBool,lBool,cList),Bg_est(bBool,lBool,cList),...
                Yif_isnan(bBool,lBool,cList),Xc]...
            = sabcondc_v3l1_gpu_batch(logYif(:,:,cList),WA(:,cList),Alibs,...
                      logT_extrap(:,:,cList),...
                      BP(:,:,cList),'lambda_a',lambda_a,'precision',precision,...
                      'Aicelib',Aicelibs,'nIter',nIter,...
                      'verbose_lad',verbose_lad,'debug_lad',debug_lad,...
                      'verbose_huwacb',verbose_huwacb,'debug_huwacb',debug_huwacb);
            Yif_cor_ori(bBool,lBool,cList) = logYif(:,:,cList) - gather(pagefun(@mtimes,gpuArray(T_est(bBool,1,cList)),gpuArray(Xc(1,:,:))));
            if ni==1
                X = nan(NA+1,nLall,nCall);
            end
            X(:,lBool,cList) = Xc;
            toc;
        end
        Yif_cor = permute(Yif_cor,[2,3,1]);
        Yif_cor_ori = permute(Yif_cor_ori,[2,3,1]);
        Yif_isnan = permute(Yif_isnan,[2,3,1]);
        Bg_est = permute(Bg_est,[2,3,1]);
        AB_est = permute(AB_est,[2,3,1]);
        T_est  = squeeze(T_est);
        X = permute(X,[2,3,1]);
        
        Yif_cor = exp(Yif_cor); T_est = exp(T_est);
        Bg_est = exp(Bg_est); AB_est = exp(AB_est);
        Yif_cor_ori = exp(Yif_cor_ori);
        
end

tend = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
fprintf('finish procceing, current time is %s.\n',tend);
fprintf('Processing time is %s\n',tend-tstart);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Now saving...\n');
% hdr, mimics CAT file production
hdr_cr = crism_const_cathdr(TRRIFdata,true);
hdr_cr.cat_history = suffix;
switch opt_img
    case 'if'
        hdr_cr.cat_input_files = [TRRIFdata.basename '.IMG'];
    case 'ra_if'
        hdr_cr.cat_input_files = [TRRIFdata.basename '.IMG'];
    case {'TRRY','TRRC','TRRB'}
        hdr_cr.cat_input_files = [basenameTRRY '.mat'];
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
    case 0
        save(fname_supple,'wa','bands','lls','T_est','ancillaries','Valid_pixels');
    case 1
        save(fname_supple,'wa','bands','lls','T_est','X','Yif_isnan');
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
    envidatawrite(single(logIce),joinPath(save_dir, [basename_Ice '.img']),hdr_cr);
    fprintf('Done\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performing interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img_sabcondl1_nan_replaced = Yif_cor;
img_modell1 = AB_est .* Bg_est;
nan_cells = isnan(Yif_cor);
img_sabcondl1_nan_replaced(nan_cells) = img_modell1(nan_cells);

Valid_pixels_good = double(Valid_pixels);
Valid_pixels_good(Valid_pixels==0) = nan;
for bi=1:nBall
    img_sabcondl1_nan_replaced(:,:,bi) = img_sabcondl1_nan_replaced(:,:,bi) .* Valid_pixels_good;
end

% replace NaN with interpolation from a model
basename_cr_nr = [basename_cr '_nr'];
hdr_cr_nr = hdr_cr;
dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
hdr_cr_nr.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced after processing.}',dt);
hdr_cr_nr.cat_history = [hdr_cr_nr.cat_history '_nr'];
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr '.hdr']));
envihdrwritex(hdr_cr_nr,joinPath(save_dir,[basename_cr_nr '.hdr']),'OPT_CMOUT',false);
fprintf('Done\n');
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr '.img']));
envidatawrite(single(img_sabcondl1_nan_replaced),joinPath(save_dir,[basename_cr_nr '.img']),hdr_cr);
fprintf('Done\n');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % gaussian filter
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logimg_sabcondl1_nan_replaced_ab = log(img_sabcondl1_nan_replaced) - log(Bg_est);
logimg_gfl = nan(size(logimg_sabcondl1_nan_replaced_ab));
gausswidth = 0.6; fltsize = 5;
h = fspecial('gaussian',fltsize,gausswidth);
for i=1:size(img_sabcondl1_nan_replaced,3)
    logimg_gfl(:,:,i) = nanimfilter(logimg_sabcondl1_nan_replaced_ab(:,:,i),h);
end
img_sabcondl1_nr_gf = exp(logimg_gfl) .* Bg_est;

basename_cr_nr_gf = [basename_cr_nr '_gf'];
hdr_cr_nr_gf = hdr_cr_nr;
dt = datetime('now','TimeZone','local','Format','eee MMM dd hh:mm:ss yyyy');
hdr_cr_nr_gf.description = sprintf('{CRISM DATA [%s] header editted timestamp, nan replaced and gauss filtered after processing modified.}',dt);
hdr_cr_nr_gf.cat_history = [hdr_cr_nr_gf.cat_history '_gf'];
hdr_cr_nr_gf.gauss_filter_std = gausswidth;
hdr_cr_nr_gf.gauss_filter_size = fltsize;
hdr_cr_nr.cat_input_files = basename_cr_nr;

fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr_gf '.hdr']));
envihdrwritex(hdr_cr_nr_gf,joinPath(save_dir,[basename_cr_nr_gf '.hdr']),'OPT_CMOUT',false);
fprintf('Done\n');
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr_gf '.img']));
envidatawrite(single(img_sabcondl1_nr_gf),joinPath(save_dir,[basename_cr_nr_gf '.img']),hdr_cr_nr_gf);
fprintf('Done\n');

fprintf('Process completed!\n');

end

