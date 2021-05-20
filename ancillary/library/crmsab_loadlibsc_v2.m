function [ Alib,infoAall,valid_idx,Aall ] = crmsab_loadlibsc_v2(optLibs,wabasename,optInterpid,c,bands_opt,wvc,cntRmvl,varargin)
global crism_env_vars
localCRISM_PDSrootDir = crism_env_vars.localCRISM_PDSrootDir;
%%
overwrite = 0;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'OVERWRITE'
                overwrite = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

optCRISMspclib = optLibs(1);
optRELAB = optLibs(2);
optUSGSsplib = optLibs(3);
optCRISMTypeLib = optLibs(4);
optInterp = crmsab_const_liboptInterp(optInterpid);
libprefix = crmsab_const_libprefix(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib);
bands = genBands(bands_opt);


Alibdir = joinPath(localCRISM_PDSrootDir,'cache/WA/',wabasename);
Alibcachefname = crmsab_constAlibcachefname(libprefix,wabasename,optInterpid,bands_opt,c);
Alibcachefilepath = joinPath(Alibdir,Alibcachefname);
% Alibcachefilepath = joinPath(Alibdir,sprintf('Alib%s_%s_interp%s_b%s_c%03d.mat',libprefix,wabasename,optInterpid,bands_opt,c));
if ~overwrite && exist(Alibcachefilepath,'file')
    if cntRmvl
        load(Alibcachefilepath,'logAcntrmvd');
    else
        load(Alibcachefilepath,'logAallNrmedvalid');
    end
    if nargout>1
        load(Alibcachefilepath,'infoAall','valid_idx');
        if nargout>3
            load(Alibcachefilepath,'Aall');
        end
    end
else
    % crism spectral library
    masterbase = crmsab_const_libmasterbase('CRISMspclib',optCRISMspclib,wabasename,optInterp(1).method,optInterp(1).retainRatio);
    cachefilepath = crmsab_const_libcachefilepath(Alibdir,masterbase,c);
    if ~exist(cachefilepath,'file')
        error('%s does not exist.\nPlease perform crism_libConvoluter first.',cachefilepath);
    else
        load(cachefilepath,'Acrismspclib','infoAcrismspclib');
        if ~isempty(Acrismspclib),Acrismspclib = Acrismspclib(bands,:);end
    end

    % RELAB
    masterbase = crmsab_const_libmasterbase('RELAB',optRELAB,wabasename,optInterp(2).method,optInterp(2).retainRatio);
    [cachefilepath] = crmsab_const_libcachefilepath(Alibdir,masterbase,c);
    if ~exist(cachefilepath,'file')
        % error('%s does not exist.\nPlease perform crism_libConvoluter first.',cachefilepath);
        fprintf('Perform crism_libConvoluter first\n');
        crism_libConvoluter('RELAB',optRELAB,'WAbasename',wabasename,'CList',c,...
            'METHOD',optInterp(2).method,'RETAINRATIO',optInterp(2).retainRatio);
    end
    load(cachefilepath,'Arelab','infoArelab');
    if ~isempty(Arelab), Arelab = Arelab(bands,:); end

    % USGS spectral library
    masterbase = crmsab_const_libmasterbase('USGSsplib',optUSGSsplib,wabasename,optInterp(3).method,optInterp(3).retainRatio);
    [cachefilepath] = crmsab_const_libcachefilepath(Alibdir,masterbase,c);
    if ~exist(cachefilepath,'file')
        error('%s does not exist.\nPlease perform crism_libConvoluter first.',cachefilepath);
    else
        load(cachefilepath,'Ausgs','infoAusgs');
        if ~isempty(Ausgs),Ausgs = Ausgs(bands,:); end
    end

    % CRISM mica library (noisy)
    masterbase = crmsab_const_libmasterbase('CRISMTypeLib',optCRISMTypeLib,wabasename,optInterp(4).method,optInterp(4).retainRatio);
    [cachefilepath] = crmsab_const_libcachefilepath(Alibdir,masterbase,c);
    if ~exist(cachefilepath,'file')
        error('%s does not exist.\nPlease perform crism_libConvoluter first.',cachefilepath);
    else
        load(cachefilepath,'AcrismTypeLib','infoAcrismTypeLib');
        if ~isempty(AcrismTypeLib), AcrismTypeLib = AcrismTypeLib(bands,:); end
    end

    Aall = [Acrismspclib Arelab Ausgs AcrismTypeLib];
    if isempty(infoAusgs)
        infoAall = merge_struct(infoAcrismspclib,infoArelab,...
                                infoAcrismTypeLib);
    else
        infoAall = merge_struct(infoAcrismspclib,infoArelab,...
                                infoAusgs,infoAcrismTypeLib);
    end
    
    % post processing:
    % prun if there is any invalid one
    valid_idx1 = ~any(isnan(Aall),1);
    valid_idx2 = all(Aall>1e-10,1);
    valid_idx = and(valid_idx1,valid_idx2);
    infoAallValid = infoAall(valid_idx);
    logAall = log(Aall);
    logAallNrmed = normalizevec(logAall,1);
    logAallNrmedvalid = logAallNrmed(:,valid_idx);
    
    [ continua,bases ] = learnContinuum(wvc,logAallNrmedvalid,...
                                        'convhull');
    [ logAcntrmvd ] = CntRmvl(logAallNrmedvalid,'additive','CONTINUA',...
                              continua );
    logAcntrmvd = -logAcntrmvd;
    
    % save
    save(Alibcachefilepath,'Aall','valid_idx','logAallNrmedvalid','logAcntrmvd',...
        'infoAall','optInterpid','optLibs','bands_opt');

end

if cntRmvl
    Alib = logAcntrmvd;
else
    Alib = logAallNrmedvalid;
end

end

