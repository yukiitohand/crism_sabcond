function [ Alib,infoAall,valid_idx,Aall_b ] = crmsab_loadAlibconv(optLibs,wabasename,optInterpid,c,bands_opt,wvc,cntRmvl,varargin)
global crism_env_vars

dir_cache = crism_env_vars.dir_CACHE;
overwrite = 0;
upd_infoAall = 0;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DIR_CACHE'
                dir_cache = varargin{i+1};
            case 'OVERWRITE'
                overwrite = varargin{i+1};
            case 'UPDATE_INFOAALL'
                upd_infoAall = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

optCRISMspclib  = optLibs(1);
optRELAB        = optLibs(2);
optUSGSsplib    = optLibs(3);
optCRISMTypeLib = optLibs(4);
libprefix = crmsab_const_libprefix(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib);

[wa_identfr] = crmsab_get_libWAIdentifier(wabasename);

Alibdir = joinPath(dir_cache,'WA',wa_identfr);
Alibcachefname = crmsab_constAlibcachefname(libprefix,wa_identfr,optInterpid,bands_opt,c);
Alibcachefilepath = joinPath(Alibdir,Alibcachefname);
if ~overwrite && exist(Alibcachefilepath,'file')
    if cntRmvl
        load(Alibcachefilepath,'logAcntrmvd');
    else
        load(Alibcachefilepath,'logAallNrmedvalid');
    end
    if nargout>1
        load(Alibcachefilepath,'valid_idx');
        infoAlibcachefname    = crmsab_constAlibcachefname_infoAall(libprefix);
        infoAlibcachefilepath = joinPath(dir_cache,infoAlibcachefname);
        load(infoAlibcachefilepath,'infoAall');
        if nargout>3
            Aall_cachefname = crmsab_constAlibcachefname_Aall(libprefix,wa_identfr,optInterpid,c);
            Aall_cachefilepath = joinPath(Alibdir,Aall_cachefname);
            load(Aall_cachefilepath,'Aall');
            propWA = crism_getProp_basenameCDR4(wabasename);
            bands = crmsab_genBands_v2(propWA.wavelength_filter,bands_opt,propWA.binning,propWA.sclk);
            Aall_b = Aall(bands,:);
        end
    end
else
    propWA = crism_getProp_basenameCDR4(wabasename);
    optInterp = crmsab_const_liboptInterp(optInterpid);
    % bands     = crmsab_genBands(bands_opt);
    bands = crmsab_genBands_v2(propWA.wavelength_filter,bands_opt,propWA.binning,propWA.sclk);
    % crism spectral library
    [Acrismspclib,infoAcrismspclib] = crmsab_load_lib_base('CRISMspclib',optCRISMspclib,wa_identfr,...
             c,'METHOD',optInterp(1).method,'retainRatio',optInterp(1).retainRatio);
    % if ~isempty(Acrismspclib), Acrismspclib = Acrismspclib(bands,:);end

    % RELAB
    [Arelab,infoArelab] = crmsab_load_lib_base('RELAB',optRELAB,wa_identfr,...
             c,'METHOD',optInterp(2).method,'retainRatio',optInterp(2).retainRatio);
    % if ~isempty(Arelab), Arelab = Arelab(bands,:); end

    % USGS spectral library
    [Ausgs,infoAusgs] = crmsab_load_lib_base('USGSsplib',optUSGSsplib,wa_identfr,...
             c,'METHOD',optInterp(3).method,'retainRatio',optInterp(3).retainRatio);
    % if ~isempty(Ausgs),Ausgs = Ausgs(bands,:); end 

    % CRISM mica library (noisy)
    [AcrismTypeLib,infoAcrismTypeLib] = crmsab_load_lib_base('CRISMTypeLib',optCRISMTypeLib,wa_identfr,...
             c,'METHOD',optInterp(4).method,'retainRatio',optInterp(4).retainRatio);
    % if ~isempty(AcrismTypeLib), AcrismTypeLib = AcrismTypeLib(bands,:); end
 
    Aall = [Acrismspclib Arelab Ausgs AcrismTypeLib];

    if ~isempty(Aall), Aall_b = Aall(bands,:); else, Aall_b = []; end

    infoAall = merge_struct(infoAcrismspclib,...
                            infoArelab,infoAusgs,infoAcrismTypeLib);
    
    % post processing:
    % prun if there is any invalid one
    valid_idx1 = ~any(isnan(Aall_b),1);
    valid_idx2 = all(Aall_b>1e-10,1);
    valid_idx = and(valid_idx1,valid_idx2);
    % infoAallValid = infoAall(valid_idx);
    logAall = log(Aall_b);
    logAallNrmed = normalizevec(logAall,1);
    logAallNrmedvalid = logAallNrmed(:,valid_idx);
    
    [ continua,bases ] = learnContinuum(wvc,logAallNrmedvalid,...
                                        'convhull');
    [ logAcntrmvd ] = CntRmvl(logAallNrmedvalid,'additive','CONTINUA',...
                              continua );
    logAcntrmvd = -logAcntrmvd;
    
    % save
    save(Alibcachefilepath,'valid_idx','logAallNrmedvalid','logAcntrmvd');
    
    % save infoAall
    infoAlibcachefname    = crmsab_constAlibcachefname_infoAall(libprefix);
    infoAlibcachefilepath = joinPath(dir_cache,infoAlibcachefname);
    if ~exist(infoAlibcachefilepath,'file') ||  upd_infoAall
        save(infoAlibcachefilepath,'infoAall');
    end
    
    % save Aall
    Aall_cachefname    = crmsab_constAlibcachefname_Aall(libprefix,wa_identfr,optInterpid,c);
    Aall_cachefilepath = joinPath(Alibdir,Aall_cachefname);
    save(Aall_cachefilepath,'Aall');

end

if cntRmvl
    Alib = logAcntrmvd;
else
    Alib = logAallNrmedvalid;
end

end

