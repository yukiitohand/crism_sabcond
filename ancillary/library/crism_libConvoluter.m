function [] = crism_libConvoluter(libname,opt,varargin)
% this is a main function to perform convolution of the different kinds of
% CRISM libraries with all of the specified types of WA files. 
% The output file is stored in 
%     pdir_cache = joinPath(localCRISM_PDSrootDir, 'cache/WA/')
% filename is defined by
%     pdir_cache2 = joinPath(pdir_cache,wabasename);
%     [masterbase] = crmsab_const_libmasterbase(libname,opt,wabasename,method,retainRatio);
%     [cachefilepath] = crmsab_const_libcachefilepath(pdir_cache2,masterbase,c);
% See those files for details
%   
%   Input Parameters
%      libname: {'CRISMspclib','RELAB','USGSsplib','CRISMTypeLib',
%                'abscoeffH2Oicelib_Mastrapa2009'
%                'abscoeffH2Oicelib_Grundy1998'
%                'abscoeffH2Oicelib_GhoSSTGrundy1998'
%                'abscoeffH2Oicelib_Warren2008'
%                'abscoeffCO2icelib_Hansen'
%                'absxsecH2Olib_HITRAN'
%                'absxsecCO2lib_HITRAN'
%                'absxsecCOlib_HITRAN'} 
%               case sensitive
%      opt: index that specify subclass of 'libname'. Please see below
%           for what subclass is defined for each 'libname'
%           'CRISMspclib' : load_CRISMspclib.m
%           'RELAB' : load_spclib_relab.m
%           'USGSsplib' : load_splibUSGS.m
%           'CRISMTypeLib' : load_CRISMTypeLib.m
%           'abscoeffH2Oicelib_Mastrapa2009': load_abscoeffH2Oicelib_Mastrapa2009.m
%           'abscoeffH2Oicelib_Grundy1998'  : load_abscoeffH2Oicelib_Grundy1998.m
%           'abscoeffH2Oicelib_GhoSSTGrundy1998': load_abscoeffH2Oicelib_GhoSSTGrundy1998.m
%           'abscoeffH2Oicelib_Warren2008'  : load_abscoeffH2Oicelib_Warren2008.m
%           'abscoeffCO2icelib_Hansen'      : load_abscoeffCO2icelib_Hansen.m
%           'absxsecH2Olib_HITRAN'          : load_absxsecH2Olib_HITRAN.m
%           'absxsecCO2lib_HITRAN'          : load_absxsecCO2lib_HITRAN.m
%           'absxsecCOlib_HITRAN'          : load_absxsecCOlib_HITRAN.m
%   Optional Parameters
%      'METHOD' : {'interpCRISMspc', 'interp1'}
%                 (default) 'interpCRISMspc'
%                 interpCRISMspc is recommended for high spectral resolutio
%                 data, otherwise use 'interp1'
%      'RETAINRATIO': option for 'interpCRISMspc.m'
%                     (default) 0.1
%      'BINNING' : pattern for regular expression, binning mode of WA file.
%                 20th character of the basename of WA {'0','1','2','3'}
%                 0: 1x, 1: 2x, 2: 5x, 3: 10x
%                 (default) '0'
%      'VERSION'  : pattern for regular expression, version of the product
%                 (default) '3'
%      'SENSOR_ID' : pattern for regular expression, sensor_id
%                    (default) 'L'
%      'OVERWRITE' : whether or not to overwrite the file
%                    (default) 0
%      'WARN_OVERWRITE': whether or not you want to get prompt whe you come
%                       across overwriting
%                       (default) 1
%      'WAbasename': string or cell array, provide if you want to perform convolution on a
%                    specific WA frame. If this is specified, then the
%                    other option of the pattern is not active.
%                    (default) ''
%      'CList'     : list of the columns to operate
%                    if it's empty, then all the columns are performed
%                    (default) []
%   Output
%      none
%   crism_libConvoluter('CRISMTypeLib',2,'METHOD','interp1','WV_BIN','3');
global crism_env_vars
dir_cache = crism_env_vars.dir_CACHE;

if nargin<2
    fprintf('Usage: [] = crism_libConvoluter(libname,opt,varargin)\n');
    return;
end
    
overwrite = 0; 
warn_overwrite = 1;
retainRatio = 0.1;
method = 'interpCRISMspc';
binning = '0';
sensor_id = 'L';
vr = '3';
wabasename = '';
cList = [];

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DIR_CACHE'
                dir_cache = varargin{i+1};
            case 'OVERWRITE'
                overwrite = varargin{i+1};
            case 'WARN_OVERWRITE'
                warn_overwrite = varargin{i+1};
            case 'METHOD'
                method = varargin{i+1};
            case 'RETAINRATIO'
                retainRatio = varargin{i+1};
            case 'BINNING'
                binning = varargin{i+1};
            case 'SENSOR_ID'
                sensor_id = varargin{i+1};
            case 'VERSION'
                vr = varargin{i+1};
            case 'WABASENAME'
                wabasename = varargin{i+1};
            case 'CLIST'
                cList = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

switch method
    case 'interpCRISMspc'
        methodbl = 1;
    case 'interp1'
        methodbl = 0;
    otherwise
        error('method %s is not valid. Choose "interpCRISMspc" or "interp1" (case sensitive).',method);
end

% read library
switch libname
    case 'CRISMspclib'
        % [CRISMspclib] = readCRISMspclib();
        % libs_CRISMspclib = 'all';
        [CRISMspclib,libs_CRISMspclib] = load_CRISMspclib(opt);
        xmult = 1;
    case 'RELAB'
        [spclib_relab_actv] = load_spclib_relab(opt);
        xmult = 1;
    case 'USGSsplib'
        [splibUSGS] = load_splibUSGS(opt);
        xmult = 1000;
    case 'CRISMTypeLib'
        [crismTypeLib] = load_CRISMTypeLib(opt);
        xmult = 1000;
    case 'abscoeffH2Oicelib_Mastrapa2009'
        [abscoeffH2Oicelib_Mastrapa2009] = load_abscoeffH2Oicelib_Mastrapa2009();
        xmult = 1;
    case 'abscoeffH2Oicelib_Grundy1998'
        [abscoeffH2Oicelib_Grundy1998] = load_abscoeffH2Oicelib_Grundy1998();
        xmult = 1;
    case 'abscoeffH2Oicelib_GhoSSTGrundy1998'
        [abscoeffH2Oicelib_GhoSSTGrundy1998] = load_abscoeffH2Oicelib_GhoSSTGrundy1998();
        xmult = 1;
    case 'abscoeffH2Oicelib_Warren2008'
        [abscoeffH2Oicelib_Warren2008] = load_abscoeffH2Oicelib_Warren2008();
        xmult = 1;
    case 'abscoeffCO2icelib_Hansen'
        [abscoeffCO2icelib_Hansen] = load_abscoeffCO2icelib_Hansen();
        xmult = 1;
    case 'absxsecH2Olib_HITRAN'
        [absxsecH2Olib_HITRAN] = load_absxsecH2Olib_HITRAN(opt);
        xmult = 1000;
    case 'absxsecCO2lib_HITRAN'
        [absxsecCO2lib_HITRAN] = load_absxsecCO2lib_HITRAN(opt);
        xmult = 1000;
    case 'absxsecCOlib_HITRAN'
        [absxsecCOlib_HITRAN] = load_absxsecCOlib_HITRAN(opt);
        xmult = 1000;
    case 'Q_H2Oicelib_Grundy1998'
        [Q_H2Oicelib_Grundy1998] = load_Q_H2Oice_Grundy1998(opt);
        xmult = 1;
    otherwise
        error('lib %s is not defined. Please select from\n{"CRISMspclib","RELAB","USGSsplib","CRISMTypLib"} (case sensitive)',libname);
end

% setups.retainRatio = retainRatio;
dir_cacheWA = joinPath(dir_cache, 'WA/');



if isempty(wabasename)
    [ propWAptr ] = crism_create_propCDR4basename( 'Acro','WA','BINNING',binning,'SENSOR_ID',sensor_id,'Version',vr);
    [~,wabasename,~] = crism_search_cdr_fromProp(propWAptr);
    if isempty(wabasename)
        error('No matching WA file is found');
    end
end


if ischar(wabasename)
    WAbasenameList = {wabasename};
elseif iscell(wabasename)
    % Eliminate duplicated WA files.
    WAbasenameList = []; wa_identfr_list = [];
    for i=1:length(wabasename)
        wabasename_i = wabasename{i};
        wa_identfr_i = crmsab_get_libWAIdentifier(wabasename_i);
        if ~any(strcmpi(wa_identfr_i,wa_identfr_list))
            WAbasenameList  = [WAbasenameList {wabasename_i}];
            wa_identfr_list = [wa_identfr_list {wa_identfr_i}];
        end
    end
end



for i=1:length(WAbasenameList)
    wabasename = WAbasenameList{i};
    [wa_identfr] = crmsab_get_libWAIdentifier(wabasename);
    %propWA = getProp_basenameCDR4(wabasename);
    %WAdir = get_dirpath_cdr_fromProp(propWA);
    WAdata = CRISMdata(wabasename,'');
    % [lblwa,hdrwa,imgwa] = crismCDRread_v2(WAdir,'WA',wabasename);
    propSB = WAdata.prop; propSB.acro_calibration_type = 'SB';
    % SBdir = get_dirpath_cdr_fromProp(propSB); 
    sbbasename = get_basenameCDR4_fromProp(propSB);
    SBdata = CRISMdata(sbbasename,'');
    % [lblsb,hdrsb,imgsb] = crismCDRread_v2(SBdir,'SB',sbbasename);
    
    imgwa = WAdata.readimgi();
    imgsb = SBdata.readimgi();
    
    dir_cacheWAbase = joinPath(dir_cacheWA,wa_identfr);
    if ~exist(dir_cacheWAbase,'dir')
        mkdir(dir_cacheWAbase);
    end
    [masterbase] = crmsab_const_libmasterbase(libname,opt,wa_identfr,method,retainRatio);
    fprintf('Starting %s, current time: %s\n',wabasename,datetime());
    if isempty(cList)
        cList = 1:WAdata.hdr.samples;
    elseif i>1
        cList = 1:WAdata.hdr.samples;
    end
    
    for c=cList
        wvc = squeeze(imgwa(:,c,:));
        if ~all(isnan(wvc))
            tc = tic;
            sbc = squeeze(imgsb(:,c,:))';
            [cachefilepath] = crmsab_const_libcachefilepath(dir_cacheWAbase,masterbase,c);
            
            if ~overwrite && exist(cachefilepath,'file')
                fprintf('Skipping %s\n',cachefilepath);
            elseif ~exist(cachefilepath,'file') || overwrite
                flg = 1;
                if exist(cachefilepath,'file')
                    if warn_overwrite
                        flg = input(sprintf('Do you want to overwrite %s?(1[yes]/0[no])',cachefilepath));
                        if flg
                            fprintf('Overwriting %s\n',cachefilepath);
                        end
                    end
                end
                if flg
                    switch libname
                        case 'CRISMspclib'
                            [Acrismspclib,infoAcrismspclib] = convCRISMspclib_v2(CRISMspclib,...
                                libs_CRISMspclib,wvc,sbc,'RETAINRATIO',retainRatio,'XMULT',xmult);
                            option.method = method; option.retainRatio = retainRatio;
                            save(cachefilepath,'Acrismspclib','infoAcrismspclib','option');
                        case 'RELAB'
                            [Arelab,option] = libstruct_convoluter(spclib_relab_actv,wvc,...
                                                    methodbl,'SB',sbc,'retainRatio',retainRatio,'XMULT',xmult);
                            infoArelab = spclib_relab_actv;
                            save(cachefilepath,'Arelab','infoArelab','option');
                        case 'USGSsplib'
                            [Ausgs,option] = libstruct_convoluter(splibUSGS,wvc,methodbl,'XMULT',xmult);
                            infoAusgs = splibUSGS;
                            save(cachefilepath,'Ausgs','infoAusgs','option');
                        case 'CRISMTypeLib'
                            [AcrismTypeLib,option] = libstruct_convoluter(crismTypeLib,wvc,methodbl,'XMULT',xmult);
                            infoAcrismTypeLib = crismTypeLib;
                            save(cachefilepath,'AcrismTypeLib','infoAcrismTypeLib','option');
                        case 'abscoeffH2Oicelib_Mastrapa2009'
                            [AabscoeffMastrapa2009,option]...
                                = libstruct_convoluter(abscoeffH2Oicelib_Mastrapa2009,wvc,...
                                methodbl,'YFieldName','abscoeff','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult);
                            infoMastrapa2009 = abscoeffH2Oicelib_Mastrapa2009;
                            save(cachefilepath,'AabscoeffMastrapa2009','infoMastrapa2009','option');
                        case 'abscoeffH2Oicelib_Grundy1998'
                            [AabscoeffGrundy1998,option]...
                                = libstruct_convoluter(abscoeffH2Oicelib_Grundy1998,wvc,...
                                methodbl,'YFieldName','abscoeff','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult);
                            infoGrundy1998 = abscoeffH2Oicelib_Grundy1998;
                            save(cachefilepath,'AabscoeffGrundy1998','infoGrundy1998','option');
                        case 'abscoeffH2Oicelib_GhoSSTGrundy1998'
                            [AabscoeffGhoSSTGrundy1998,option]...
                                = libstruct_convoluter(abscoeffH2Oicelib_GhoSSTGrundy1998,wvc,...
                                methodbl,'YFieldName','abscoeff','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult);
                            infoGhoSSTGrundy1998 = abscoeffH2Oicelib_GhoSSTGrundy1998;
                            save(cachefilepath,'AabscoeffGhoSSTGrundy1998','infoGhoSSTGrundy1998','option');
                        case 'abscoeffH2Oicelib_Warren2008'
                            [AabscoeffWarren2008,option]...
                                = libstruct_convoluter(abscoeffH2Oicelib_Warren2008,wvc,...
                                methodbl,'YFieldName','abscoeff','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult);
                            infoWarren2008 = abscoeffH2Oicelib_Warren2008;
                            save(cachefilepath,'AabscoeffWarren2008','infoWarren2008','option');
                        case 'abscoeffCO2icelib_Hansen'
                            [AabscoeffHansen,option]...
                                = libstruct_convoluter(abscoeffCO2icelib_Hansen,wvc,...
                                methodbl,'YFieldName','abscoeff','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult);
                            infoHansen = abscoeffCO2icelib_Hansen;
                            save(cachefilepath,'AabscoeffHansen','infoHansen','option');
                        case 'absxsecH2Olib_HITRAN'
                            [AabsxsecH2Olib_HITRAN,option]...
                                = libstruct_convoluter(absxsecH2Olib_HITRAN,wvc,...
                                methodbl,'YFieldName','xsec','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult,'Batch',1);
                            infoabsxsecH2Olib_HITRAN = absxsecH2Olib_HITRAN;
                            save(cachefilepath,'AabsxsecH2Olib_HITRAN','infoabsxsecH2Olib_HITRAN','option');
                        case 'absxsecCO2lib_HITRAN'
                            [AabsxsecCO2lib_HITRAN,option]...
                                = libstruct_convoluter(absxsecCO2lib_HITRAN,wvc,...
                                methodbl,'YFieldName','xsec','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult,'Batch',1);
                            infoabsxsecCO2lib_HITRAN = absxsecCO2lib_HITRAN;
                            save(cachefilepath,'AabsxsecCO2lib_HITRAN','infoabsxsecCO2lib_HITRAN','option');
                        case 'absxsecCOlib_HITRAN'
                            [AabsxsecCOlib_HITRAN,option]...
                                = libstruct_convoluter(absxsecCOlib_HITRAN,wvc,...
                                methodbl,'YFieldName','xsec','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult,'Batch',1);
                            infoabsxsecCOlib_HITRAN = absxsecCOlib_HITRAN;
                            save(cachefilepath,'AabsxsecCOlib_HITRAN','infoabsxsecCOlib_HITRAN','option');
                        case 'Q_H2Oicelib_Grundy1998'
                            [AQ_H2Oicelib_Grundy1998,option]...
                                = libstruct_convoluter(Q_H2Oicelib_Grundy1998,wvc,...
                                methodbl,'YFieldName','Q','SB',sbc,'retainRatio',retainRatio,'XMULT',xmult,'Batch',1);
                            infoQ_H2Oicelib_Grundy1998 = Q_H2Oicelib_Grundy1998;
                            save(cachefilepath,'AQ_H2Oicelib_Grundy1998','infoQ_H2Oicelib_Grundy1998','option');
                        otherwise
                            error('lib %s is not defined',libname);
                    end
                    
                end
            end
            ts = toc(tc);
            fprintf('column %03d finished %f[s]\n',c,ts);
        end
    end
    fprintf('Finshed %s, current time: %s\n',wabasename, datetime());
    
end

end

