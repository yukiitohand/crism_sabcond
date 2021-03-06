function [Yif_cor,Bg_est,AB_est,T_est] = sabcondv3_Te(obs_id,varargin)
% huwacb option
nIter = 5;
vis = 0;
lambda_a = 0.0045;
% band option
opt_img = 'TRRY';
bands_opt = 4;
% library options
optCRISMspclib = 1;
optRELAB = 1;
optUSGSsplib = 1;
optCRISMTypeLib = 2;
% opticelib = 1;
% opthitranlib = 1;
cntRmvl = 1;
optInterpid = 1;
t_mode = 2;
% file name option
mt = 'sabcondv3'; % method type
additional_suffix = 'l1_gadmm_a_v2_Te';
isdebug = false;
save_pdir = './resu/';
save_dir_yyyy_doy = false;
force = false;
gausssigma = 0.6;
optBP = 'self'; %{'pri','all','none','self'}
skip_ifexist = false;
Lib_Modify = '0';

global localCRISM_PDSrootDir

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'NITER'
                nIter = varargin{i+1};
            case 'VIS'
                vis = varargin{i+1};
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
            %case 'OPT_ICELIB'
            %    opticelib = varargin{i+1};
            %case 'OPT_HITRANLIB'
            %    opthitranlib = varargin{i+1};
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
            case 'DEBUG'
                isdebug = varargin{i+1};
            case 'LIB_MODIFY'
                Lib_Modify = varargin{i+1};
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
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
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
%fprintf('t_mode:%d\n',t_mode);
fprintf('gauss_sigma:%f\n',gausssigma);
fprintf('opt_img: %s\n',opt_img);
fprintf('optBP: %s\n',optBP);

bands = crmsab_genBands(bands_opt);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,'','');
% libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,opthitranlib);

% open log file
username = char(java.lang.System.getProperty('user.name'));
fname = sprintf('log_%s_%s.txt',username,datetime('now','TimeZone','local','Format','yyyyMMdd'));
diary(joinPath(save_pdir,fname));

fprintf('Current directory:%s\n',pwd);

%% Read image and ancillary data and format them for processing
crism_obs = CRISMObservationFRT(obs_id,'SENSOR_ID','L');
if ~isempty(crism_obs.info.basenameIF)
    crism_obs.load_data(crism_obs.info.basenameIF,crism_obs.info.dir_trdr,'if');
    TRRIFdata = crism_obs.data.if;
    TRRIF_is_empty = false;
else
    TRRIFdata = '';
    TRRIF_is_empty = true;
end
if ~isempty(crism_obs.info.basenameRA)
    crism_obs.load_data(crism_obs.info.basenameRA,crism_obs.info.dir_trdr,'ra');
    TRRRAdata = crism_obs.data.ra;
else
    TRRRAdata = '';
end
if ~isempty(crism_obs.info.basenameSC)
    crism_obs.load_data(crism_obs.info.basenameSC,crism_obs.info.dir_edr,'sc');
    EDRdata = crism_obs.data.sc;
else
    EDRdata = '';
end

if TRRIF_is_empty
    TRRIFdata = TRRRAdata;
end

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


propIF = getProp_basenameOBSERVATION(TRRIFdata.basename);
if TRRIF_is_empty
    propIF.activity_id = 'IF';
end
basenameIF = get_basenameOBS_fromProp(propIF);

basename_cr = [basenameIF suffix];
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
SBdata = TRRIFdata.readCDR('SB'); SBdata.readimgi();
% crim = CRISMImage(obs_id,'SENSOR_ID','L');
nLall = TRRIFdata.hdr.lines; nCall = TRRIFdata.hdr.samples; nBall = TRRIFdata.hdr.bands;
% lines = 2:nLall-1;
llines = 1:nLall;
nL = length(llines);
nB = length(bands);

lBool = false(nLall,1); lBool(llines) = true;
bBool = false(nBall,1); bBool(bands) = true;


basenameWA = WAdata.basename;
% SB = permute(SBdata.img(:,:,bands),[3,1,2]);
WA = squeeze(WAdata.img(:,:,bands))';

switch opt_img
    case 'if'
        if TRRIF_is_empty
            error('TRR I/F does not exist.');
        else
            Yif = TRRIFdata.readimgi();
        end
    case 'ra_if'
        TRRRAIFdata = crism_obs.load_data([crism_obs.info.basenameRA '_IF'],crism_obs.info.dir_trdr,'ra_if');
        Yif = TRRRAIFdata.readimgi();
    case 'yuki_IoF'
        d_IoF = joinPath(localCRISM_PDSrootDir,'./../YUK/', crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.product_type = 'YUK';
        prop.version = 5;
        basenameYUK2 = get_basenameOBS_fromProp(prop);
        load(joinPath(d_IoF,[basenameYUK2 '.mat']),'IoF_woc');
        Yif = IoF_woc;
        clear IoF_woc;
        Yif = flip(Yif,3);
%     case 'TRRY'
%         d_IoF = joinPath(localCRISM_PDSrootDir,'./../YUK/', crism_obs.info.yyyy_doy, crism_obs.info.dirname);
%         prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
%         prop.version = 'Y';
%         basenameTRRY = get_basenameOBS_fromProp(prop);
%         load(joinPath(d_IoF,[basenameTRRY '.mat']),'IoF_woc');
%         Yif = IoF_woc;
%         clear IoF_woc;
%         Yif = flip(Yif,3);
    case 'TRRY'
        d_IoF = joinPath(localCRISM_PDSrootDir,'./../YUK/', crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = 'Y';
        basenameTRRY = get_basenameOBS_fromProp(prop);
        TRRYIFdata = CRISMdata(basenameTRRY,d_IoF);
        Yif = TRRYIFdata.readimgi();
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

Yif = Yif(llines,:,bands);
Yif(Yif<=1e-8) = nan;
logYif = log(Yif);
logYif = permute(logYif,[3,1,2]);

clear Yif;
fprintf('finish loading Image\n');

%%
% read bad pixel
switch EDRdata.lbl.OBSERVATION_TYPE
    case {'FRT','HRL','HRS'}
        crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
        DFdata1 = crism_obs.data.df1;
        crism_obs.load_data(crism_obs.info.basenameDF{2},crism_obs.info.dir_edr,'df2');
        DFdata2 = crism_obs.data.df2;
        if strcmpi(optBP,'self')
            DFdata1.readimgi();
            DFdata2.readimgi();
            dfimage = cat(1,DFdata1.img,DFdata2.img);
            Noutliers = 4;
        end
    case {'FRS','ATO'}
        if ischar(crism_obs.info.basenameDF)
            crism_obs.load_data(crism_obs.info.basenameDF,crism_obs.info.dir_edr,'df1');
        elseif iscell(crism_obs.info.basenameDF)
            crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
        end
        if strcmpi(optBP,'self')
            DFdata1.readimgi();
            dfimage = DFdata1.img;
            Noutliers = 2;
        end
    otherwise
        error('Please define for other cases')
end

TRRIFdata.load_basenamesCDR();
% read bad pixel data
TRRIFdata.readCDR('BP');
switch EDRdata.lbl.OBSERVATION_TYPE
    case {'FRT','HRL','HRS'}
        for i=1:length(TRRIFdata.cdr.BP)
            bpdata = TRRIFdata.cdr.BP(i);
            if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                    BPdata1 = bpdata;
                elseif any(strcmpi(DFdata2.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                    BPdata2 = bpdata;
                end
            else
                BPdata_post = bpdata;
            end
        end
    case {'FRS','ATO'}
        % in case of FRS, DFdata1 and DFdata2 are same.
        for i=1:length(TRRIFdata.cdr.BP)
            bpdata = TRRIFdata.cdr.BP(i);
            if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                    BPdata1 = bpdata; BPdata2 = bpdata;
                end
            else
                BPdata_post = bpdata;
            end
        end
    otherwise
        error('Undefined observation type %s.',EDRdata.lbl.OBSERVATION_TYPE);
end



BPdata1.readimgi(); BPdata2.readimgi();
BP_pri_bool = or(BPdata1.img,BPdata2.img);
BP_pri_bool = BP_pri_bool(:,:,bands);
% BP_pri = double(BP_pri_bool);
% BP_pri(BP_pri==0) = nan;

% BP_pri = permute(BP_pri,[3,1,2]);

BPdata_post.readimgi();
% BP_all = double(BPdata_post.img(:,:,bands));
% BP_all(BP_all==0) = nan;

% BP_onlypost = double(and(BP_all==1,BP_pri_bool==0));
% BP_onlypost(BP_onlypost==0) = nan;
% BP_onlypost = permute(BP_onlypost,[3,1,2]);

GP_pri = double(BP_pri_bool==0);
GP_pri(GP_pri==0) = nan;
GP_pri = permute(GP_pri,[3,1,2]);

GP_all = double(BPdata_post.img(:,:,bands)==0);
GP_all(GP_all==0) = nan;
GP_all = permute(GP_all,[3,1,2]);

GP_self = nan(size(GP_all));
DFstd_self = nan(size(GP_all));
for c = 1:nCall
    dfimagec = squeeze(dfimage(:,c,bands));
    std_df = robust('std',dfimagec,1,'NOutliers',Noutliers);
    GPc = double(std_df<4.5);
    GPc(GPc==0) = nan;
    GP_self(:,:,c) = GPc(:);
    DFstd_self(:,:,c) = std_df(:);
end

switch optBP
    case 'pri'
        GP = GP_pri;
    case 'all'
        GP = GP_all;
    case 'self'
        GP = GP_self;
    case 'none'
        GP = true(size(GP_all));
    otherwise 
        error('optBP=%s is not defined.',optBP);
end


%% Library customization
[Alib,infoAall,valid_idx] = crmsab_loadlibsc_v2(optLibs,basenameWA,optInterpid,floor(nCall/2),bands_opt,WA(:,floor(nCall/2)),cntRmvl);
switch Lib_Modify
    case '0'
        mag_factors = ones(1,length(infoAall));
        exclude_lib_idx_all = [];
    case 'a'
        [kieserite,k_i] = searchby_multfield({'name','spc_name'},{'kieserite','crism_typespec_mono_hyd_sulf'},infoAall);
        [pyrox,pyrox_i] = searchby_multfield({'name','spc_name'},{'pyroxene'},infoAall);
        exclude_lib_idx_all = k_i;
        mag_factors = ones(1,length(infoAall));
        mag_factors(:,k_i) = 5;
        mag_factors(:,pyrox_i) = 2;
        
    case 'b'
        [kieserite,k_i] = searchby_multfield({'name','spc_name'},{'kieserite','crism_typespec_mono_hyd_sulf'},infoAall);
        [water, w_i] = searchby_multfield({'name','spc_name'},{'h2o'},infoAall);
        [gypsum, g_i] = searchby_multfield({'name','spc_name'},{'gypsum'},infoAall);
        [poly_sulf, poly_sulf_i] = searchby_multfield({'name','spc_name'},{'POLY_HYD_SULF','magnesium sulfate'},infoAall);
        [pyrox,pyrox_i] = searchby_multfield({'name','spc_name'},{'pyroxene'},infoAall);
        exclude_lib_idx_all = [k_i w_i g_i poly_sulf_i];
        mag_factors = ones(1,length(infoAall));
        mag_factors(:,k_i) = 5;
        mag_factors(:,pyrox_i) = 2;
    case 'c'
        [kieserite,k_i] = searchby_multfield({'name','spc_name'},{'kieserite','crism_typespec_mono_hyd_sulf'},infoAall);
        [pyrox,pyrox_i] = searchby_multfield({'name','spc_name'},{'pyroxene'},infoAall);
        exclude_lib_idx_all = [];
        mag_factors = ones(1,length(infoAall));
        mag_factors(:,k_i) = 5;
        mag_factors(:,pyrox_i) = 2;
    otherwise
        error('Lib_Modify %s is not defined',Lib_Modify);
end

exclude_lib_idx_all_bool = false(1,length(infoAall));
exclude_lib_idx_all_bool(exclude_lib_idx_all) = true;
        
%%
if isdebug
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % original competitors
    phot = 0; 
    atmt_src = 'trial'; %{'tbench','auto','user','default'}
    bandset_id = 'mcg'; %{'mcg','pel'}
    enable_artifact = 1;
    acro_catatp = sprintf('phot%d_%s_%s_a%d',phot,atmt_src,bandset_id,enable_artifact);
    suffix = ['_corr_' acro_catatp];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TRRIFdata_corr = crism_obs.load_data([TRRIFdata.basename suffix],TRRIFdata.dirpath,acro_catatp);
    TRRIFdata_corr.readimgi();
    Yif_cat = TRRIFdata_corr.img;

    Yif_cat = Yif_cat(llines,:,bands);
    Yif_cat(Yif_cat<=1e-8) = nan;
    logYifc_cat = log(Yif_cat);
    logYifc_cat = permute(logYifc_cat,[3,1,2]);

    clear Yif_cat
    TRRRAdata = crism_obs.data.ra;
    TRRRAIFdata_corr = crism_obs.load_data([TRRRAdata.basename '_IF' suffix],TRRRAdata.dirpath,acro_catatp);
    TRRRAIFdata_corr.readimgi();
    Yraif_cat = TRRRAIFdata_corr.img;

    Yraif_cat = Yraif_cat(llines,:,bands);
    Yraif_cat(Yraif_cat<=1e-8) = nan;
    logYraifc_cat = log(Yraif_cat);
    logYraifc_cat = permute(logYraifc_cat,[3,1,2]);
    clear Yraif_cat

%     acro_catatp_raif = ['ra_' acro_catatp];
%     suffix_raif = ['_if_corr_' acro_catatp];
%     crim.readlblhdr([crim.basename.ra suffix_raif],crim.dirpath.ra,acro_catatp_raif);
%     crim.readimg(acro_catatp_raif);
% 
%     Yraif_cat = crim.img.(acro_catatp_raif).all;
%     Yraif_cat = Yraif_cat(lines,:,bands);
%     Yraif_cat(Yraif_cat<=1e-8) = nan;
%     logYraifc_cat = log(Yraif_cat);
%     logYraifc_cat = permute(logYraifc_cat,[3,1,2]);
% 
%     clear Yraif_cat
end

%% read ADR transmission data
prop = crism_getProp_basenameCDR4(WAdata.basename);
% [ at_trans ] = load_adr( 'WV_BIN',crim.info.cdr.WA(20),'T_MODE',t_mode );
switch t_mode
    case {1,2,3}
        [ at_trans ] = crism_load_ADR_VS('t_mode',t_mode,'BINNING',prop.binning,...
                                   'WAVELENGTH_FILTER',prop.wavelength_filter);
    case {4}
        [ at_trans ] = crmsab_load_T();
    otherwise
        error('Undefined t_mode %d',t_mode);
end

T = at_trans(:,:,bands);
T(T<=1e-8) = nan;
logT = log(T);
logT = permute(logT,[3,1,2]);

clear T

fprintf('finish loading ADR\n');
%% main loop
fprintf('Start processing\n');
tstart = datetime('now','TimeZone','America/New_York','Format','d-MMM-y HH:mm:ss Z');
fprintf('Current time is %s.\n',tstart);

Yif_cor = nan([nLall,nCall,nBall]);
Yif_isnan = nan([nLall,nCall,nBall]);
T_est = nan([nBall,nCall]);
Bg_est = nan([nLall,nCall,nBall]); AB_est = nan([nLall,nCall,nBall]);
ancillaries = struct('X',cell(1,nCall),'lambda',cell(1,nCall),'nIter',cell(1,nCall),...
    'huwacb_func',cell(1,nCall),'maxiter_huwacb',cell(1,nCall),'tol_huwacb',cell(1,nCall),'gp_bool',cell(1,nCall));
RR_ori = nan([nLall,nCall,nBall]);
Yif_cor_ori = nan([nLall,nCall,nBall]);
Valid_pixels = false([nLall,nCall]);
% for c=32:nCall-9
for c = 1:nCall
%for c=290:380
% for c=[1:32,(nCall-9):ncAll]
    % ADR data is filtered, so the spectra at edges are removed.
    if ~all(isnan(WA(:,c)))
        if any(isnan(logT(:,:,c)))
            if ~any(isnan(logT(:,:,c+1)))
                logtc = logT(:,:,c+1);
            elseif ~any(isnan(logT(:,:,c-1)))
                logtc = logT(:,:,c-1);
            elseif ~any(isnan(logT(:,:,c+2)))
                logtc = logT(:,:,c+2);
            elseif ~any(isnan(logT(:,:,c-2)))
                logtc = logT(:,:,c-2);
            else
                logtc = logT(:,:,c);
            end
        else
            logtc = logT(:,:,c);
        end
        tic;
        % load library
        [Alib,infoAall,valid_idx] = crmsab_loadlibsc_v2(optLibs,basenameWA,optInterpid,c,bands_opt,WA(:,c),cntRmvl);
        mag_factors_Alib = mag_factors(valid_idx);
        exclude_lib_idx_bool_Alib = exclude_lib_idx_all_bool(valid_idx);
        Alib2 = Alib.*mag_factors_Alib;

        [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori_c,vldpxl_c]...
            = sabcondc_v3l1_e(Alib2,logYif(:,:,c),WA(:,c),logtc,'GP',GP(:,:,c),...
                'LAMBDA_A',lambda_a,'NITER',nIter,'VIS',vis,'Exclude_lib_idx',exclude_lib_idx_bool_Alib);
                %'Debug',0,'LOGYIFC_CAT',logYifc_cat(:,:,c),'LOGYRAIFC_CAT',logYraifc_cat(:,:,c));
        
        Yif_cor(lBool,c,bBool) = reshape(logYifc_cor',[nL,1,nB]);
        Yif_cor_ori(lBool,c,bBool) = reshape(logYifc_cor_ori',[nL,1,nB]);
        T_est(bBool,c) = logt_est;
        Bg_est(lBool,c,bBool) = reshape(logBg',[nL,1,nB]);
        AB_est(lBool,c,bBool) = reshape(logAB',[nL,1,nB]);
        Yif_isnan(lBool,c,bBool) = reshape(logYifc_isnan',[nL,1,nB]);
        RR_ori(lBool,c,bBool) = reshape(rr_ori_c',[nL,1,nB]);
        ancillaries(c) = ancillary;
        Valid_pixels(lBool,c) = vldpxl_c';
%         t = toc(tc);
        toc;
%         fprintf('%03d fin. %f [s]\n',c,t);
    end
end

Yif_cor = exp(Yif_cor); T_est = exp(T_est);
Bg_est = exp(Bg_est); AB_est = exp(AB_est);
Yif_cor_ori = exp(Yif_cor_ori);
GP_ori = nan([1,nCall,nBall]);
GP = permute(GP,[2,3,1]);
GP_ori(:,:,bBool) = GP;

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
    case 'yuki_IoF'
        hdr_cr.cat_input_files = [basenameYUK2 '.mat'];
    case 'TRRY'
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

fname_supple = joinPath(save_dir,[TRRIFdata.basename suffix '.mat']);
wa = WAdata.img;
wa = squeeze(wa)';
fprintf('Saving %s ...\n',fname_supple);
save(fname_supple,'wa','bands','llines','T_est','ancillaries','Valid_pixels','GP_ori');
fprintf('Done\n');

basename_Bg = [TRRIFdata.basename suffix '_Bg'];
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg '.hdr']));
envihdrwritex(hdr_cr,joinPath(save_dir,[basename_Bg '.hdr']),'OPT_CMOUT',false);
fprintf('Done\n');
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_Bg '.img']));
envidatawrite(single(Bg_est),joinPath(save_dir, [basename_Bg '.img']),hdr_cr);
fprintf('Done\n');


basename_AB = [TRRIFdata.basename suffix '_AB'];
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB '.hdr']));
envihdrwritex(hdr_cr,joinPath(save_dir, [basename_AB '.hdr']),'OPT_CMOUT',false);
fprintf('Done\n');
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_AB '.img']));
envidatawrite(single(AB_est),joinPath(save_dir, [basename_AB '.img']),hdr_cr);
fprintf('Done\n');

% residual
basename_RR = [TRRIFdata.basename suffix '_RR'];
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_RR '.hdr']));
envihdrwritex(hdr_cr,joinPath(save_dir, [basename_RR '.hdr']),'OPT_CMOUT',false);
fprintf('Done\n');
fprintf('Saving %s ...\n',joinPath(save_dir, [basename_RR '.img']));
envidatawrite(single(RR_ori),joinPath(save_dir, [basename_RR '.img']),hdr_cr);
fprintf('Done\n');

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

% logimg_sabcondl1_nan_replaced_ab = log(img_sabcondl1_nan_replaced) - log(Bg_est);
% logimg_gfl = nan([nLall,nCall,nBall]);
% for i=1:nBall
%     logimg_gfl(:,:,i) = imgaussfilt(logimg_sabcondl1_nan_replaced_ab(:,:,i),gausssigma);
% end
% img_sabcondl1_nr_gf = exp(logimg_gfl + log(Bg_est));
% 
% basename_cr_nr_gf = [basename_cr_nr '_gf'];
% fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr_gf '.hdr']));
% envihdrwritex(hdr_cr,joinPath(save_dir,[basename_cr_nr_gf '.hdr']),'OPT_CMOUT',false);
% fprintf('Done\n');
% fprintf('Saving %s ...\n',joinPath(save_dir, [basename_cr_nr_gf '.img']));
% envidatawrite(single(img_sabcondl1_nr_gf),joinPath(save_dir,[basename_cr_nr_gf '.img']),hdr_cr);
% fprintf('Done\n');

fprintf('Process completed!\n');

end

