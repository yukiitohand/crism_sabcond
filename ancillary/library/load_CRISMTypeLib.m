function [crismTypeLib] = load_CRISMTypeLib(opt_CRISMTypeLib,varargin)
% load MRO CRISM Type Spectra Library
% Input Parameters
%   opt_CRISMTypeLib : subclass id
% Output Parameters
%   crismTypeLib: struct, see 'readCRISMTypeLibrary.m' for detail
%         additional field
%            'spclib': 'CRISM Type Spectral Library'
%            'cumidex'
%            'reflectannce'

global crism_env_vars
dir_cache = crism_env_vars.dir_CACHE;

overwrite = 0; 
warn_overwrite = 1;
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
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

[infocachefname] = sprintf('CRISMTypeLib%d_info.mat',opt_CRISMTypeLib);
infocachefilepath = fullfile(dir_cache,infocachefname);
if ~overwrite && exist(infocachefilepath,'file')
    load(infocachefilepath,'infoAcrismTypeLib');
    crismTypeLib = infoAcrismTypeLib;
elseif ~exist(infocachefilepath,'file') || overwrite
    flg=1;
    if exist(infocachefilepath,'file')
        if warn_overwrite
            flg = doyouwantto('overwrite',sprintf('%s exists',infocachefilepath));
            if flg
                fprintf('Overwriting %s\n',infocachefilepath);
            end
        end
    end
    
    if flg
        switch opt_CRISMTypeLib
            case 1
                % scale and shift the ratioed reflectances.
                [crismTypeLib] = readCRISMTypeLibrary();
                load('crismTypeLib_scale.mat','sList','b1List','b2List',...
                     'refList');
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = ...
                    (crismTypeLib(i).ratioed_cor-b1List(i))*sList(i)+b2List(i);
                end
            case 2
                % the original ratioed reflectances are used
                [crismTypeLib] = readCRISMTypeLibrary();
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            case {0,3}
                % formaly case 3 (changed on Dec. 11, 2018)
                % empty
                [crismTypeLib] = [];
            case {4}
                % (added on Jan. 16, 2019)
                % exclude ices.
                [crismTypeLib] = readCRISMTypeLibrary();
                [~,ices_i] = searchby_multfield('name','ice',crismTypeLib);
                nonices_i = setdiff(1:length(crismTypeLib),ices_i);
                crismTypelib_woIce = crismTypeLib(nonices_i);
                crismTypeLib = crismTypelib_woIce;
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            case {5}
                % (added on June 8, 2022)
                % include only h2o_ice
                [crismTypeLib] = readCRISMTypeLibrary();
                [~,h2oices_i] = searchby_multfield('name','h2o_ice',crismTypeLib);
                crismTypelib_H2OIce = crismTypeLib(h2oices_i);
                crismTypeLib = crismTypelib_H2OIce;
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            case {6}
                % (added on Aug. 23, 2022)
                % exclude gypsum
                [crismTypeLib] = readCRISMTypeLibrary();
                [~,gypsum_i] = searchby_multfield('name','crism_typespec_gypsum',crismTypeLib);
                nongypsum_i = setdiff(1:length(crismTypeLib),gypsum_i);
                crismTypeLib = crismTypeLib(nongypsum_i);
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            case {7}
                % (added on Aug. 23, 2022)
                % exclude gypsum and ices
                [crismTypeLib] = readCRISMTypeLibrary();
                [~,gypsum_i] = searchby_multfield('name','crism_typespec_gypsum',crismTypeLib);
                [~,ices_i] = searchby_multfield('name','ice',crismTypeLib);
                nongypice_i = setdiff(1:length(crismTypeLib),[gypsum_i ices_i]);
                crismTypeLib = crismTypeLib(nongypice_i);
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            case {8}
                % (added on Sep. 2, 2022)
                % include only h2o_ice and co2_ice
                [crismTypeLib] = readCRISMTypeLibrary();
                [~,ices_i] = searchby_multfield('name','ice',crismTypeLib);
                crismTypelib_Ice = crismTypeLib(ices_i);
                crismTypeLib = crismTypelib_Ice;
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            case {9}
                % (added on Sep. 2, 2022)
                % include only co2_ice
                [crismTypeLib] = readCRISMTypeLibrary();
                [~,co2ices_i] = searchby_multfield('name','co2_ice',crismTypeLib);
                crismTypelib_CO2Ice = crismTypeLib(co2ices_i);
                crismTypeLib = crismTypelib_CO2Ice;
                for i=1:length(crismTypeLib)
                    crismTypeLib(i).reflectance = crismTypeLib(i).ratioed_cor;
                end
            otherwise
                error('opt_CRISMTypeLib %d is not defined',opt_CRISMTypeLib);
        end
        
        infoAcrismTypeLib = crismTypeLib;
        save(infocachefilepath,'infoAcrismTypeLib');
    end
end

if ~isempty(crismTypeLib)
    [crismTypeLib.spclib] = deal('CRISM Type Spectral Library');
    iList = num2cell(1:length(crismTypeLib));
    [crismTypeLib.cumindex] = iList{:};
end

end