function [Yif_cor,Bg_est,AB_est,T_est,Yif_cor_nr,Yif_cor_ori,...
    Yif_isnan,RR_ori,Valid_pixels,ancillary,infoAall,valid_idx,GP]...
    = sabcondv3_e_columntest(obs_id,cList,varargin)
% huwacb option
nIter = 5;
vis = 0;
lambda_a = 0.01;
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
isdebug = false;
gausssigma = 0.6;
optBP = 'pri'; %{'pri','all','none'}
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
            case 'LIB_MODIFY'
                Lib_Modify = varargin{i+1};
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
            case 'OPTINTERPID'
                optInterpid = varargin{i+1};
            case 'DEBUG'
                isdebug = varargin{i+1};
            case 'OPTBP'
                optBP = varargin{i+1};

            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


fprintf('lambda_a:%f\n',lambda_a);
fprintf('nIter:%d\n',nIter);
fprintf('t_mode:%d\n',t_mode);
fprintf('gauss_sigma:%f\n',gausssigma);
fprintf('opt_img: %s\n',opt_img);
fprintf('optBP: %s\n',optBP);

bands = crmsab_genBands(bands_opt);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];

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

%% Read image and ancillary data and format them for processing
TRRIFdata.load_basenamesCDR();
WAdata = TRRIFdata.readCDR('WA'); WAdata.readimgi();
SBdata = TRRIFdata.readCDR('SB'); SBdata.readimgi();
% crim = CRISMImage(obs_id,'SENSOR_ID','L');
nLall = TRRIFdata.hdr.lines; nCall = TRRIFdata.hdr.samples; nBall = TRRIFdata.hdr.bands;
% lines = 2:nLall-1;
lls = 1:nLall;
nL = length(lls);
nB = length(bands);

lBool = false(nLall,1); lBool(lls) = true;
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
            %Yif = TRRIFdata.lazyEnviReadci(c);
        end
    case 'ra_if'
        TRRRAIFdata = crism_obs.load_data([crism_obs.info.basenameRA '_IF'],crism_obs.info.dir_trdr,'ra_if');
        %Yif = TRRRAIFdata.lazyEnviReadci(c);
        Yif = TRRRAIFdata.readimgi();
    case 'TRRY'
        d_IoF = joinPath(localCRISM_PDSrootDir,'./../YUK/', crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = 'Y';
        basenameTRRY = get_basenameOBS_fromProp(prop);
        TRRYIFdata = CRISMdata(basenameTRRY,d_IoF);
        Yif = TRRYIFdata.readimgi();
        %Yif = TRRYIFdata.lazyEnviReadci(c);
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

Yif = Yif(lls,cList,bands);
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
    case {'FRS','ATO'}
        if ischar(crism_obs.info.basenameDF)
            crism_obs.load_data(crism_obs.info.basenameDF,crism_obs.info.dir_edr,'df1');
        elseif iscell(crism_obs.info.basenameDF)
            crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
        end
        DFdata1 = crism_obs.data.df1;
        DFdata2 = crism_obs.data.df1;
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

switch optBP
    case 'pri'
        GP = GP_pri;
    case 'all'
        GP = GP_all;
    case 'none'
        GP = true(size(GP_all));
    otherwise 
        error('optBP=%s is not defined.',optBP);
end

%% Library customization
[Alib,infoAall,valid_idx] = crmsab_loadlibsc_v2(optLibs,basenameWA,optInterpid,floor(cList(1)),bands_opt,WA(:,floor(cList(1))),cntRmvl);
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
    case 'd'
        [kieserite,k_i] = searchby_multfield({'name','spc_name'},{'kieserite','crism_typespec_mono_hyd_sulf'},infoAall);
        [water, w_i] = searchby_multfield({'name','spc_name'},{'h2o'},infoAall);
        [gypsum, g_i] = searchby_multfield({'name','spc_name'},{'gypsum'},infoAall);
        [poly_sulf, poly_sulf_i] = searchby_multfield({'name','spc_name'},{'POLY_HYD_SULF','magnesium sulfate'},infoAall);
        [szomolnokite,szo_i] = searchby_multfield({'name','spc_name'},'szomolnokite',infoAall);
        [pyrox,pyrox_i] = searchby_multfield({'name','spc_name'},{'pyroxene'},infoAall);
        exclude_lib_idx_all = [k_i w_i g_i poly_sulf_i szo_i];
        mag_factors = ones(1,length(infoAall));
        mag_factors(:,k_i) = 5;
        mag_factors(:,pyrox_i) = 2;
        mag_factors(:,szo_i) = 5;
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

    Yif_cat = Yif_cat(lls,:,bands);
    Yif_cat(Yif_cat<=1e-8) = nan;
    logYifc_cat = log(Yif_cat);
    logYifc_cat = permute(logYifc_cat,[3,1,2]);

    clear Yif_cat
    TRRRAdata = crism_obs.data.ra;
    TRRRAIFdata_corr = crism_obs.load_data([TRRRAdata.basename '_IF' suffix],TRRRAdata.dirpath,acro_catatp);
    TRRRAIFdata_corr.readimgi();
    Yraif_cat = TRRRAIFdata_corr.img;

    Yraif_cat = Yraif_cat(lls,:,bands);
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
nCLength = length(cList);
Yif_cor = nan([nBall,nLall,nCLength]);
Yif_cor_nr = nan([nBall,nLall,nCLength]);
Yif_isnan = nan([nBall,nLall,nCLength]);
T_est = nan([nBall,nCLength]);
Bg_est = nan([nBall,nLall,nCLength]); AB_est = nan([nBall,nLall,nCLength]);
RR_ori = nan([nBall,nLall,nCLength]);
Yif_cor_ori = nan([nBall,nLall,nCLength]);
Valid_pixels = false([nLall,nCLength]);
for ci=1:nCLength
    c = cList(ci);
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

        [Alib,infoAall,valid_idx] = crmsab_loadlibsc_v2(optLibs,basenameWA,optInterpid,c,bands_opt,WA(:,c),cntRmvl);
        mag_factors_Alib = mag_factors(valid_idx);
        exclude_lib_idx_bool_Alib = exclude_lib_idx_all_bool(valid_idx);
        Alib2 = Alib.*mag_factors_Alib;

        [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori_c,vldpxl_c]...
            = sabcondc_v3l1_e(Alib2,logYif(:,:,ci),WA(:,c),logtc,'GP',GP(:,:,c),...
                'LAMBDA_A',lambda_a,'NITER',nIter,'VIS',vis,'Exclude_lib_idx',exclude_lib_idx_bool_Alib);
    %                 'LOGYIFC_CAT',logYifc_cat(:,:,c));%'LOGYRAIFC_CAT',logYraifc_cat(:,:,c));
        logmodel = logBg + logAB;
        logYif_cor_nr = logYifc_cor;
        logYif_cor_nr(logYifc_isnan) = logmodel(logYifc_isnan);

        Yif_cor(bBool,lBool,ci) = logYifc_cor;
        Yif_cor_nr(bBool,lBool,ci) = logYif_cor_nr;
        Yif_cor_ori(bBool,lBool,ci) = logYifc_cor_ori;
        T_est(bBool,ci) = logt_est;
        Bg_est(bBool,lBool,ci) = logBg;
        AB_est(bBool,lBool,ci) = logAB;
        Yif_isnan(bBool,lBool,ci) = logYifc_isnan;
        RR_ori(bBool,lBool,ci) = rr_ori_c;
        Valid_pixels(lBool,ci) = vldpxl_c;

        toc;
    else
        error('Maybe c=% 3d is out of range',c);
    end
end

Yif_cor = exp(Yif_cor); T_est = exp(T_est);
Bg_est = exp(Bg_est); AB_est = exp(AB_est);
Yif_cor_ori = exp(Yif_cor_ori);
Yif_cor_nr = exp(Yif_cor_nr);
GP = GP(:,:,cList);


end

