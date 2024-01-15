function [Aicelib,infoAicelib,valid_idx] = crmsab_loadAicelibconv(opt,wabasename,c,bands_opt,wvc,varargin)

global crism_env_vars

dir_cache = crism_env_vars.dir_CACHE;
overwrite = 0;
cntRmvl   = 0;
upd_infoAicelib = 0;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DIR_CACHE'
                dir_cache = varargin{i+1};
            case 'CNTRMVL'
                cntRmvl = varargin{i+1};
            case 'OVERWRITE'
                overwrite = varargin{i+1};
            case 'UPDATE_INFOAAICELIB'
                upd_infoAicelib = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

[wa_identfr] = crmsab_get_libWAIdentifier(wabasename);

Aicelibdir = joinPath(dir_cache,'WA/',wa_identfr);
Aicelibcachefname = crmsab_constAicelibcachefname(wa_identfr,opt,bands_opt,c,cntRmvl);
Aicelibcachefilepath = joinPath(Aicelibdir,Aicelibcachefname);
if ~overwrite && exist(Aicelibcachefilepath,'file')
    load(Aicelibcachefilepath,'Aicelib');
    if nargout>1
        load(Aicelibcachefilepath,'valid_idx');
        infoAicelibcachefname    = crmsab_constAicelibcachefname_infoAicelib(opt);
        infoAicelibcachefilepath = joinPath(dir_cache,infoAicelibcachefname);
        load(infoAicelibcachefilepath,'infoAicelib');
    end
else
    propWA = crism_getProp_basenameCDR4(wabasename);
    [Aicelib,infoAicelib] = crmsab_load_icelibconv_wrapper(opt,wa_identfr,c);
    % select the band option
    bands = crmsab_genBands_v2(propWA.wavelength_filter,bands_opt,propWA.binning,propWA.sclk);
    Aicelib= Aicelib(bands,:); 
    % prun if there is any invalid one
    valid_idx = ~any(isnan(Aicelib),1);
%     valid_idx2 = all(Aicelib>1e-10,1);
%     valid_idx = and(valid_idx1,valid_idx2);
    infoAicelib = infoAicelib(valid_idx);
    Aicelib = Aicelib(:,valid_idx);
    if cntRmvl
        [ continua,bases ] = learnContinuum(wvc,Aicelib,'convhull');
        [ Aicelib ] = CntRmvl(Aicelib,'additive','CONTINUA',continua );
        Aicelib = -Aicelib;
    end
    
    save(Aicelibcachefilepath,'Aicelib','valid_idx');
    infoAicelibcachefname    = crmsab_constAicelibcachefname_infoAicelib(opt);
    infoAicelibcachefilepath = joinPath(dir_cache,infoAicelibcachefname);
    if ~exist(infoAicelibcachefilepath,'file') || upd_infoAicelib
        save(infoAicelibcachefilepath,'infoAicelib');
    end
end


end

