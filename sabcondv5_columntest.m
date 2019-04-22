function [Yif_cor,Bg_est,AB_est,T_est,Yif_cor_nr,Yif_cor_ori,...
    Yif_isnan,RR_ori,Valid_pixels,ancillary,infoAall,valid_idx,GP]...
    = sabcondv5_columntest(obs_id,cList,varargin)
% huwacb option
nIter = 5;
vis = 0;
lambda_a = 0.01;
% band option
opt_img = 'TRRB';
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
% file name option
isdebug = false;
gausssigma = 0.6;
optBP = 'none'; %{'pri','all','none'}
lls = [];
T_add = [];

global crism_env_vars
dir_yuk = crism_env_vars.dir_YUK;

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
            case 'T_MODE'
                t_mode = varargin{i+1};
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
            case 'OPTINTERPID'
                optInterpid = varargin{i+1};
            case 'LLS'
                lls = varargin{i+1};
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
fprintf('gauss_sigma:%f\n',gausssigma);
fprintf('opt_img: %s\n',opt_img);
fprintf('optBP: %s\n',optBP);

bands = genBands(bands_opt);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
fprintf('Current directory:%s\n',pwd);

%% load crism_observation and TRR I/F data
crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');

TRRIFdata = load_CRISMdata(crism_obs.info.basenameIF,crism_obs.info.dir_trdr);
TRRRAdata = load_CRISMdata(crism_obs.info.basenameRA,crism_obs.info.dir_trdr);
EDRdata = load_CRISMdata(crism_obs.info.basenameSC,crism_obs.info.dir_edr);

if isempty(TRRIFdata)
    TRRIFdata = TRRRAdata;
end
if isempty(TRRIFdata)
    error('TRRIF data seems to be missing.');
end

%% Read image and format it for processing
% read some attributed information: wavelength information, band, 
% line information is computed.
TRRIFdata.load_basenamesCDR();
WAdata = TRRIFdata.readCDR('WA'); WAdata.readimgi();
SBdata = TRRIFdata.readCDR('SB'); SBdata.readimgi();
nLall = TRRIFdata.hdr.lines; nCall = TRRIFdata.hdr.samples; nBall = TRRIFdata.hdr.bands;
if isempty(lls)
    lls = 1:nLall;
end
nL = length(lls); nB = length(bands);
lBool = false(nLall,1); lBool(lls) = true;
bBool = false(nBall,1); bBool(bands) = true;

basenameWA = WAdata.basename;
WA = squeeze(WAdata.img(:,:,bands))';

d_yuk_trr = joinPath(dir_yuk,crism_obs.info.yyyy_doy,crism_obs.info.dirname);
switch lower(opt_img)
    case 'if'
        if isempty(TRRIFdata)
            error('TRR I/F does not exist.');
        else
            Yif = TRRIFdata.readimgi();
        end
    case 'ra_if'
        TRRRAIFdata = CRISMdata([crism_obs.info.basenameRA '_IF'],crism_obs.info.dir_trdr);
        Yif = TRRRAIFdata.readimgi();
    case 'trry'
        trr_ver = 'Y';
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = trr_ver;
        basenameTRRY = get_basenameOBS_fromProp(prop);
        TRRYIFdata = CRISMdata(basenameTRRY,d_yuk_trr);
        Yif = TRRYIFdata.readimgi();
    case 'trrb'
        trr_vr = 'B';
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = trr_vr;
        basenameTRRB = get_basenameOBS_fromProp(prop);
        TRRBIFdata = CRISMdata(basenameTRRB,d_yuk_trr);
        Yif = TRRBIFdata.readimgi();
    case 'trrc'
        trr_vr = 'C';
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = trr_vr;
        basenameTRRB = get_basenameOBS_fromProp(prop);
        TRRBIFdata = CRISMdata(basenameTRRB,d_yuk_trr);
        Yif = TRRBIFdata.readimgi();
    otherwise
        error('opt_img = %s is not defined',opt_img);
end

Yif = Yif(lls,:,bands);
Yif(Yif<=1e-8) = nan;
Yif = permute(Yif,[3,1,2]);
logYif = log(Yif);

fprintf('finish loading Image\n');

%% Bad pixel information
TRRIFdata.load_basenamesCDR();
TRRIFdata.readCDR('BP');
[BPdata1,BPdata2,BPdata_post] = load_BPdata(TRRIFdata);

[BP_pri1nan] = formatBPpri1nan(BPdata1,BPdata2,'band_inverse',true);
[GP_pri] = convertBP1nan2GP1nan(BP_pri1nan);
GP_pri = permute(GP_pri,[1,3,2]);
GP_pri = GP_pri(bands,:,:);
[BP_post1nan] = formatBP1nan(BPdata_post,'band_inverse',true);
[GP_all] = convertBP1nan2GP1nan(BP_post1nan);
GP_all = permute(GP_all,[1,3,2]);
GP_all = GP_all(bands,:,:);

BP_none = all(isnan(Yif),2);
GP_none = double(~BP_none);
GP_none(GP_none==0) = nan;
% GP_none = permute(GP_none,[1,3,2]);

switch optBP
    case 'pri'
        GP = GP_pri;
    case 'all'
        GP = GP_all;
    case 'self'
        GP = GP_self;
    case 'none'
        GP = GP_none;
    otherwise 
        error('optBP=%s is not defined.',optBP);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the weight for each dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch EDRdata.lbl.OBSERVATION_TYPE
    case {'FRT','HRL','HRS'}
        DFdata1 = CRISMdata(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr);
        DFdata2 = CRISMdata(crism_obs.info.basenameDF{2},crism_obs.info.dir_edr);
        Noutliers = 4;
    case {'FRS','ATO'}
        if ischar(crism_obs.info.basenameDF)
            DFdata1 = CRISMdata(crism_obs.info.basenameDF,crism_obs.info.dir_edr);
        elseif iscell(crism_obs.info.basenameDF)
            DFdata1 = CRISMdata(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr);
        end
        Noutliers = 2;
        DFdata2 = [];
    otherwise
        error('Please define for other cases')
end

% load processed dark files
switch EDRdata.lbl.OBSERVATION_TYPE
    case {'FRT','HRL','HRS'}
        propDF1_IF = DFdata1.prop;
        propDF1_IF.activity_id = 'IF';
        propDF1_IF.product_type = 'TRR';
        propDF1_IF.version = trr_vr;
        bnameDF1_IF = get_basenameOBS_fromProp(propDF1_IF);
        load(joinPath(d_yuk_trr,[bnameDF1_IF '.mat']),'IoF_bk1_o');
        IoF_bk1_o = flip(IoF_bk1_o,3);
        IoF_bk1_o = permute(IoF_bk1_o,[3,1,2]);
        
        propDF2_IF = DFdata2.prop;
        propDF2_IF.activity_id = 'IF';
        propDF2_IF.product_type = 'TRR';
        propDF2_IF.version = trr_vr;
        bnameDF2_IF = get_basenameOBS_fromProp(propDF2_IF);
        load(joinPath(d_yuk_trr,[bnameDF2_IF '.mat']),'IoF_bk2_o');
        IoF_bk2_o = flip(IoF_bk2_o,3);
        IoF_bk2_o = permute(IoF_bk2_o,[3,1,2]);
        
    case {'FRS','ATO'}
        propDF1_IF = DFdata1.prop;
        propDF1_IF.activity_id = 'IF';
        propDF1_IF.product_type = 'TRR';
        propDF1_IF.version = trr_vr;
        bnameDF1_IF = get_basenameOBS_fromProp(propDF1_IF);
        load(joinPath(d_yuk_trr,[bnameDF1_IF '.mat']),'IoF_bk1_o');
        IoF_bk1_o = flip(IoF_bk1_o,3);
        IoF_bk1_o = permute(IoF_bk1_o,[3,1,2]);
    otherwise
        error('Please define for other cases')
end

ifdfstd_self = nan(size(GP_all));
ifdfImage = cat(2,IoF_bk1_o,IoF_bk2_o);
for cd = 1:nCall
    ifdfimagec = ifdfImage(bands,:,cd);
    % std_df = robust_v3('stdl1',ifdfimagec,2,'NOutliers',Noutliers);
    std_df = robust_v3('med_abs_dev_from_med',ifdfimagec,2,'NOutliers',Noutliers);
    ifdfstd_self(:,:,cd) = std_df(:);
end


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

end

%% read ADR transmission data
propWA = getProp_basenameCDR4(WAdata.basename);
% [ at_trans ] = load_adr( 'WV_BIN',crim.info.cdr.WA(20),'T_MODE',t_mode );
switch t_mode
    case {1,2,3}
        [ at_trans ] = load_ADR_VS('BINNING',propWA.binning,...
                                   'WAVELENGTH_FILTER',propWA.wavelength_filter);
    case {4}
        [ at_trans ] = load_T();
    otherwise
        error('Undefined t_mode %d',t_mode);
end

T = at_trans(:,:,bands);
T(T<=1e-8) = nan;
T = permute(T,[3,1,2]);
logT = log(T);
% clear T

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

        [Alib,infoAall,valid_idx] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,bands_opt,WA(:,c),cntRmvl);

        [ logt_est,logYifc_cor,logAB,logBg,logYifc_cor_ori,logYifc_isnan,ancillary,rr_ori_c,vldpxl_c]...
            = sabcondc_v5l1_med(Alib,logYif(:,:,c),WA(:,c),logtc,'GP',GP(:,:,c),...
              'LAMBDA_A',lambda_a,'NITER',nIter,'VIS',vis,'T',T,'Yif',Yif(:,:,c),'stdl1_ifdf',ifdfstd_self(:,:,c),...
              'debug',true,'BP_pri',BP_pri1nan(bands,c),'BP_All',BP_post1nan(bands,c));
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

function [CRISMdataobj] = load_CRISMdata(basename,dirpath)
    if ~isempty(basename)
        CRISMdataobj = CRISMdata(basename,dirpath);
    else
        CRISMdataobj = [];
    end
end
