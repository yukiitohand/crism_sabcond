function [Alib,infoA,option] = crmsab_load_lib_base(libname,opt,wabasename,c,varargin)
% Load spectral librrary convolved to CRISM wavelength channels
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
%                'abscoeffH2Oicelib_Warren2008'} 
%               case sensitive
%      opt: index that specify subclass of 'libname'. Please see below
%           for what subclass is defined for each 'libname'
%           'CRISMspclib' : load_CRISMspclib.m
%           'RELAB' : load_spclib_relab.m
%           'USGSsplib' : load_splibUSGS.m
%           'CRISMTypeLib' : load_CRISMTypeLib.m
%           'abscoeffH2Oicelib_Mastrapa2009': load_abscoeffH2Oicelib_Mastrapa2009.m
%           'abscoeffH2Oicelib_Grundy1998'  : load_abscoeffH2Oicelib_Grundy1998.m
%           'abscoeffH2Oicelib_Warren2008'  : load_abscoeffH2Oicelib_Warren2008.m
%           'absxsecH2Olib_HITRAN'          
%      wabasename: basename of WA files for which the library is convolved
%      c: specify column number
%   Optional Parameters
%      'METHOD' : {'interpCRISMspc', 'interp1'}
%                 (default) 'interpCRISMspc'
%                 interpCRISMspc is recommended for high spectral resolutio
%                 data, otherwise use 'interp1'
%      'RETAINRATIO': option for 'interpCRISMspc.m'
%                     (default) 0.1
%   Output
%      Alib: library matrix [L x N], L is the number of wavelength samples,
%            N is the number of spectra in the library
%      infoA: N-length struct, storing detailed information for the library
%             elements.
%      option: option for the 
%   crism_libConvoluter('CRISMTypeLib',2,'METHOD','interp1','WV_BIN','3');

global crism_env_vars
dir_cache = crism_env_vars.dir_CACHE;

if nargin<4
    fprintf('Usage: [Alib,infoA,option] = crism_load_lib(libname,opt,wabasename,c,varargin)\n');
    return;
end
    
retainRatio = 0.1;
method = 'interpCRISMspc';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'METHOD'
                method = varargin{i+1};
            case 'RETAINRATIO'
                retainRatio = varargin{i+1};
            case 'DIR_CACHE'
                dir_cache = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

% switch method
%     case 'interpCRISMspc'
%         methodbl = 1;
%     case 'interp1'
%         methodbl = 0;
%     otherwise
%         error('method %s is not valid. Choose "interpCRISMspc" or "interp1" (case sensitive).',method);
% end

if opt
    switch libname
        case 'CRISMspclib'
            Alibname = 'Acrismspclib'; infoAname = 'infoAcrismspclib';
        case 'RELAB'
            Alibname = 'Arelab'; infoAname = 'infoArelab';
        case 'USGSsplib'
            Alibname = 'Ausgs'; infoAname = 'infoAusgs';
        case 'CRISMTypeLib'
            Alibname = 'AcrismTypeLib'; infoAname = 'infoAcrismTypeLib';
        case 'abscoeffH2Oicelib_Mastrapa2009'
            Alibname = 'AabscoeffMastrapa2009'; infoAname = 'infoMastrapa2009';
        case 'abscoeffH2Oicelib_Grundy1998'
            Alibname = 'AabscoeffGrundy1998'; infoAname = 'infoGrundy1998';
        case 'abscoeffH2Oicelib_GhoSSTGrundy1998'
            Alibname = 'AabscoeffGhoSSTGrundy1998'; infoAname = 'infoGhoSSTGrundy1998';
        case 'abscoeffH2Oicelib_Warren2008'
            Alibname = 'AabscoeffWarren2008'; infoAname = 'infoWarren2008';
        case 'abscoeffCO2icelib_Hansen'
            Alibname = 'AabscoeffHansen'; infoAname = 'infoHansen';
        case 'absxsecH2Olib_HITRAN'
            Alibname = 'AabsxsecH2Olib_HITRAN'; infoAname = 'infoabsxsecH2Olib_HITRAN';
        case 'absxsecCO2lib_HITRAN'
            Alibname = 'AabsxsecCO2lib_HITRAN'; infoAname = 'infoabsxsecCO2lib_HITRAN';
        case 'absxsecCOlib_HITRAN'
            Alibname = 'AabsxsecCOlib_HITRAN'; infoAname = 'infoabsxsecCOlib_HITRAN';
        case 'Q_H2Oicelib_Grundy1998'
            Alibname = 'AQ_H2Oicelib_Grundy1998'; infoAname = 'infoQ_H2Oicelib_Grundy1998';
        case 'Hapke_H2Oicelib_Grundy1998'
            Alibname = 'AHapke_H2Oicelib_Grundy1998'; infoAname = 'infoHapke_H2Oicelib_Grundy1998';
        otherwise
            error([ ...
                'lib %s is not defined. Please select from: \n'           , ...
                '     CRISMspclib\n'                                      , ...
                '     RELAB\n'                                            , ...
                '     USGSsplib\n'                                        , ...
                '     CRISMTypLib\n'                                      , ...
                '     abscoeffH2Oicelib_Mastrapa2009\n'                   , ...
                '     abscoeffH2Oicelib_Grundy1998\n'                     , ...
                '     abscoeffH2Oicelib_GhoSSTGrundy1998\n'               , ...
                '     abscoeffH2Oicelib_Warren2008\n'                     , ...
                '     abscoeffCO2icelib_Hansen\n'                         , ...
                '     absxsecH2Olib_HITRAN\n'                             , ...
                '     absxsecCO2lib_HITRAN\n'                             , ...
                '     absxsecCOlib_HITRAN\n'                              , ...
                '     Q_H2Oicelib_Grundy1998\n'                           , ...
                '(case sensitive)'                                          ...
                ], ...
                libname);
    end

    dir_cacheWA = joinPath(dir_cache,'WA',wabasename);
    [libcachefname]    = crmsab_const_libcachefname(libname,opt,wabasename,method,retainRatio,c);
    [libcachefilepath] = joinPath(dir_cacheWA,libcachefname);


    switch nargout
        case 1
            data = load(libcachefilepath,Alibname);
            Alib = data.(Alibname);
        case 2
            data = load(libcachefilepath,Alibname);
            Alib = data.(Alibname);
            [info_cachefname]    = crmsab_const_libcachefname_info(libname,opt);
            [info_cachefilepath] = joinPath(dir_cache,info_cachefname);
            data = load(info_cachefilepath,infoAname);
            infoA = data.(infoAname);
        case 3
            data = load(libcachefilepath,Alibname,'option');
            Alib = data.(Alibname);
            option = data.option;
            [info_cachefname]    = crmsab_const_libcachefname_info(libname,opt);
            [info_cachefilepath] = joinPath(dir_cache,info_cachefname);
            data = load(info_cachefilepath,infoAname);
            infoA = data.(infoAname);
        otherwise
            error('The number of output arguments is invalid');
    end

else
    Alib = []; infoA = []; option = [];
end

end


