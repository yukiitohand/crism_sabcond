function [ Alib,infoAall,valid_idx ] = crmsab_loadAlibconv_mini(optLibs,wabasename,optInterpid,c,bands_opt,varargin)
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

% cntRmvl = 1;

optCRISMspclib  = optLibs(1);
optRELAB        = optLibs(2);
optUSGSsplib    = optLibs(3);
optCRISMTypeLib = optLibs(4);
libprefix = crmsab_const_libprefix(optCRISMspclib,optRELAB,optUSGSsplib,optCRISMTypeLib);

[wa_identfr] = crmsab_get_libWAIdentifier(wabasename);

Alibdir = joinPath(dir_cache,'WA',wa_identfr);
Alibcachefname = crmsab_constAlibcachefname_mini(libprefix,wa_identfr,optInterpid,bands_opt,c);
Alibcachefilepath = joinPath(Alibdir,Alibcachefname);
if exist(Alibcachefilepath,'file')
    load(Alibcachefilepath,'logAcntrmvd');
    if nargout>1
        load(Alibcachefilepath,'valid_idx');
        infoAlibcachefname    = crmsab_constAlibcachefname_infoAall(libprefix);
        infoAlibcachefilepath = joinPath(dir_cache,infoAlibcachefname);
        load(infoAlibcachefilepath,'infoAall');
    end
else
    error('Library %s is found in %s.',Alibcachefname,Alibdir);
end

Alib = logAcntrmvd;

end

