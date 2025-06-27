function [ Asurficelib,infoAsurficelib,valid_idx ] = ...
    crmsab_loadAsurficelibconv(opt, wabasename, c, bands_opt, wvc, ...
    varargin)
%
% Read 
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

Asurficelibdir = joinPath(dir_cache,'WA/',wa_identfr);
Asurficelibcachefname = crmsab_constAsurficelibcachefname(wa_identfr,opt,bands_opt,c,cntRmvl);
Asurficelibcachefilepath = joinPath(Asurficelibdir,Asurficelibcachefname);
if ~overwrite && exist(Asurficelibcachefilepath,'file')
    if cntRmvl
        load(Asurficelibcachefilepath,'logAsurficelib_cntrmvd');
    else
        load(Asurficelibcachefilepath,'logAsurficelibNrmed');
    end
    if nargout>1
        load(Asurficelibcachefilepath,'valid_idx');
        infoAsurficelibcachefname    = crmsab_constAsurficelibcachefname_infoAsurficelib(opt);
        infoAsurficelibcachefilepath = joinPath(dir_cache,infoAsurficelibcachefname);
        load(infoAsurficelibcachefilepath,'infoAsurficelib');
    end
else
    propWA = crism_getProp_basenameCDR4(wabasename);
    [Asurficelib,infoAsurficelib] = crmsab_load_surficelibconv_wrapper(opt,wa_identfr,c);
    % select the band option
    bands = crmsab_genBands_v2(propWA.wavelength_filter,bands_opt,propWA.binning,propWA.sclk);
    Asurficelib = Asurficelib(bands,:); 
    % prun if there is any invalid one
    valid_idx1 = ~any(isnan(Asurficelib),1);
    valid_idx2 = all(Asurficelib>1e-10,1);
    valid_idx = and(valid_idx1,valid_idx2);
    infoAsurficelib = infoAsurficelib(valid_idx);
    Asurficelib = Asurficelib(:,valid_idx);

    logAsurficelib = log(Asurficelib);
    logAsurficelibNrmed = normalizevec(logAsurficelib,1);

    if cntRmvl
        [ continua,bases ] = learnContinuum(wvc,logAsurficelibNrmed,'convhull');
        [ logAsurficelib_cntrmvd ] = CntRmvl(logAsurficelibNrmed,'additive','CONTINUA',continua );
        logAsurficelib_cntrmvd = -logAsurficelib_cntrmvd;
    end
    
    save(Asurficelibcachefilepath,'logAsurficelib_cntrmvd','logAsurficelibNrmed', 'valid_idx');
    infoAsurficelibcachefname    = crmsab_constAsurficelibcachefname_infoAsurficelib(opt);
    infoAsurficelibcachefilepath = joinPath(dir_cache,infoAsurficelibcachefname);
    if ~exist(infoAsurficelibcachefilepath,'file') || upd_infoAicelib
        save(infoAsurficelibcachefilepath,'infoAsurficelib');
    end

end

if cntRmvl
    Asurficelib = logAsurficelib_cntrmvd;
else
    Asurficelib = logAsurficelibNrmed;
end


end
