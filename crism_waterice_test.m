function [water_ice_result,water_ice_exist,Xice_mean_insig,Xice_map,insig_map,ancillary] = crism_waterice_test(obs_id,varargin)
% huwacb option
nIter = 5;
vis = 0;
lambda_a = 0.01;
% band option
opt_img = 'TRRY';
bands_opt = 3;
% library options
optCRISMspclib = 1;
optRELAB = 1;
optUSGSsplib = 1;
optCRISMTypeLib = 2;
opticelib = 1;
cntRmvl = 1;
optInterpid = 1;
t_mode = 2;
optBP = 'pri'; %{'pri','all','none'}

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
            case 'OPT_ICELIB'
                opticelib = varargin{i+1};
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
            case 'TRRY_PDIR'
                dir_yuk = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

fprintf('lambda_a:%f\n',lambda_a);
fprintf('nIter:%d\n',nIter);
fprintf('t_mode:%d\n',t_mode);
fprintf('opt_img: %s\n',opt_img);
fprintf('optBP: %s\n',optBP);
fprintf('opt_icelib: %d\n',opticelib);

bands = genBands(bands_opt);
optLibs = [optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib];
%libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,'','');
libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,'');
% libprefix = const_libprefix_v2(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib,opticelib,opthitranlib);

% open log file
%username = char(java.lang.System.getProperty('user.name'));
%fname = sprintf('log_%s_%s.txt',username,datetime('now','TimeZone','local','Format','yyyyMMdd'));
%diary(joinPath(save_pdir,fname));

fprintf('Current directory:%s\n',pwd);

%% Read image and ancillary data and format them for processing
crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');
% crism_obsS = CRISMObservation(obs_id,'SENSOR_ID','S');
switch upper(crism_obs.info.obs_classType)
    case {'FRT','HRL','HRS','FRS','ATO','MSP','HSP'}
        if ~isempty(crism_obs.info.basenameIF)
            TRRIF_is_empty = false;
            TRRIFdata = CRISMdata(crism_obs.info.basenameIF,'');
        elseif ~isempty(crism_obs.info.basenameRA)
            TRRIF_is_empty = true;
            TRRIFdata = CRISMdata(crism_obs.info.basenameRA,'');
        else
            error('Check data');
        end
    case {'FFC'}
        switch ffc_counter
            case 1
                if ~isempty(crism_obs.info.basenameIF)
                    TRRIF_is_empty = false;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameIF{1},'');
                elseif ~isempty(crism_obs.info.basenameRA)
                    TRRIF_is_empty = true;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameRA{1},'');
                else
                    error('Check data');
                end
            case 3
                if ~isempty(crism_obs.info.basenameIF)
                    TRRIF_is_empty = false;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameIF{2},'');
                elseif ~isempty(crism_obs.info.basenameRA)
                    TRRIF_is_empty = true;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameRA{2},'');
                else
                    error('Check data');
                end
            otherwise
                error('Check data');
        end
end
[DFdata1,DFdata2] = get_DFdata4SC(TRRIFdata,crism_obs);
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
%     case 'yuki_IoF'
%         d_IoF = joinPath(dir_yuk, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
%         prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
%         prop.product_type = 'YUK';
%         prop.version = 5;
%         basenameYUK2 = get_basenameOBS_fromProp(prop);
%         load(joinPath(d_IoF,[basenameYUK2 '.mat']),'IoF_woc');
%         Yif = IoF_woc;
%         clear IoF_woc;
%         Yif = flip(Yif,3);
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
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = 'Y';
        basenameTRRY = get_basenameOBS_fromProp(prop);
        d_IoF = joinPath(dir_yuk, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        TRRYIFdata = CRISMdata(basenameTRRY,d_IoF);
        Yif = TRRYIFdata.readimgi();
    case 'TRRB'
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = 'B';
        basenameTRRB = get_basenameOBS_fromProp(prop);
        d_IoF = joinPath(dir_yuk, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        TRRBIFdata = CRISMdata(basenameTRRB,d_IoF);
        Yif = TRRBIFdata.readimgi();
    case 'TRRC'
        prop = getProp_basenameOBSERVATION(TRRIFdata.basename);
        prop.version = 'C';
        basenameTRRC= get_basenameOBS_fromProp(prop);
        d_IoF = joinPath(dir_yuk, crism_obs.info.yyyy_doy, crism_obs.info.dirname);
        TRRCIFdata = CRISMdata(basenameTRRC,d_IoF);
        Yif = TRRCIFdata.readimgi();
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
[BPdata1,BPdata2,BPdata_post] = load_BPdataSC_fromDF(TRRIFdata,DFdata1.basename,DFdata2.basename);

% switch EDRdata.lbl.OBSERVATION_TYPE
%     case {'FRT','HRL','HRS'}
%         crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
%         DFdata1 = crism_obs.data.df1;
%         crism_obs.load_data(crism_obs.info.basenameDF{2},crism_obs.info.dir_edr,'df2');
%         DFdata2 = crism_obs.data.df2;
%     case {'FRS','ATO'}
%         crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
%         DFdata1 = crism_obs.data.df1;
%         DFdata2 = crism_obs.data.df1;
%     otherwise
%         error('Please define for other cases')
% end

% TRRIFdata.load_basenamesCDR();
% EDRdata = crism_obs.data.sc;
% % read bad pixel data
% TRRIFdata.readCDR('BP');
% switch EDRdata.lbl.OBSERVATION_TYPE
%     case {'FRT','HRL','HRS'}
%         for i=1:length(TRRIFdata.cdr.BP)
%             bpdata = TRRIFdata.cdr.BP(i);
%             if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                 if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                     BPdata1 = bpdata;
%                 elseif any(strcmpi(DFdata2.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                     BPdata2 = bpdata;
%                 end
%             else
%                 BPdata_post = bpdata;
%             end
%         end
%     case {'FRS','ATO'}
%         % in case of FRS, DFdata1 and DFdata2 are same.
%         for i=1:length(TRRIFdata.cdr.BP)
%             bpdata = TRRIFdata.cdr.BP(i);
%             if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                 if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                     BPdata1 = bpdata; BPdata2 = bpdata;
%                 end
%             else
%                 BPdata_post = bpdata;
%             end
%         end
%     otherwise
%         error('Undefined observation type %s.',EDRdata.lbl.OBSERVATION_TYPE);
% end

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

%% read ADR transmission data
% prop = getProp_basenameCDR4(WAdata.basename);
% [ at_trans ] = load_adr( 'WV_BIN',crim.info.cdr.WA(20),'T_MODE',t_mode );
switch t_mode
    case {1,2,3}
        [ at_trans ] = load_ADR_VS('BINNING',WAdata.prop.binning,...
                                   'WAVELENGTH_FILTER',WAdata.prop.wavelength_filter);
    case {4}
        [ at_trans ] = load_T();
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

water_ice_exist = nan([1,nCall,1]);
Xice_mean_insig = nan([1,nCall,1]);
Xice_map = nan([nLall,nCall]);
insig_map = nan([nLall,nCall]);

ancillaries = struct('X',cell(1,nCall),'lambda',cell(1,nCall),'nIter',cell(1,nCall),...
    'huwacb_func',cell(1,nCall),'maxiter_huwacb',cell(1,nCall),'tol_huwacb',cell(1,nCall),'gp_bool',cell(1,nCall));

for c = 1:50:nCall
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

        [Alib,infoAall,valid_idx] = loadlibsc_v2(optLibs,basenameWA,optInterpid,c,bands_opt,WA(:,c),cntRmvl);
        [Aicelib,infoAicelib] = loadlibc_crism_icelib(opticelib,basenameWA,c,bands_opt,WA(:,c),'overwrite',0,'CNTRMVL',0);

        [ water_ice_existc,Xice1_mean_insig, Xice1, idx_insig,ancillary]...
            = sabcondc_v4l1_b_icetest(Alib,Aicelib,logYif(:,:,c),WA(:,c),logtc,'GP',GP(:,:,c),...
              'LAMBDA_A',lambda_a,'Lambda_a_ice',0.000,'NITER',nIter,'VIS',vis,'maxiter',200);
          
        ancillaries(c) = ancillary;
        Xice_map(lBool,c) = Xice1';
        insig_map(lBool,c) = idx_insig';
        water_ice_exist(c) = water_ice_existc;
        Xice_mean_insig(c) = Xice1_mean_insig;
        toc;
%         fprintf('%03d fin. %f [s]\n',c,t);
    end
end

if nanmean(water_ice_exist)>0.8
    water_ice_result = 1;
else
    water_ice_result = 0;
end

end

