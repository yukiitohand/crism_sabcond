function [Aicelib,infoAicelib,valid_idx] = crmsab_loadAicelibconv_mini(opt,wabasename,c,bands_opt,varargin)

global crism_env_vars

dir_cache = crism_env_vars.dir_CACHE;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DIR_CACHE'
                dir_cache = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

cntRmvl = 0;

[wa_identfr] = crmsab_get_libWAIdentifier(wabasename);

Aicelibdir = joinPath(dir_cache,'WA/',wa_identfr);
Aicelibcachefname = crmsab_constAicelibcachefname(wa_identfr,opt,bands_opt,c,cntRmvl);
Aicelibcachefilepath = joinPath(Aicelibdir,Aicelibcachefname);
if ~overwrite && exist(Aicelibcachefilepath,'file')
    load(Aicelibcachefilepath,'Aicelib');
    if nargout>1
        load(Aicelibcachefilepath,'valid_idx');
        infoAicelibcachefname    = crmsab_constAicelibcachefname_infoAicelib(libprefix);
        infoAicelibcachefilepath = joinPath(dir_cache,infoAicelibcachefname);
        load(infoAicelibcachefilepath,'infoAicelib');
    end
else
    error('Library %s is found in %s.',Aicelibcachefname,Aicelibdir);

end


end

