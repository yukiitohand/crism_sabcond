function [Hapke_H2Oicelib_Grundy1998] = load_Hapke_H2Oice_Grundy1998(opt,varargin)
% [abscoeffH2Oicelib_Grundy1998] = load_abscoeffH2Oicelib_Grundy1998()
%   load the library of Hapke's bidirectional reflectance of H2O ice computed from the
%   optical constants published in (Grundy et al., 1998).
%
% Output Parameters
%   Hapke_H2Oicelib_Grundy1998: struct with 
%   

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

[infocachefname] = sprintf('Hapke_H2Oicelib_Grundy1998%d_info.mat',opt);
infocachefilepath = fullfile(dir_cache,infocachefname);
if ~overwrite && exist(infocachefilepath,'file')
    load(infocachefilepath,'infoHapke_H2Oicelib_Grundy1998');
    Hapke_H2Oicelib_Grundy1998 = infoHapke_H2Oicelib_Grundy1998;
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

        [abscoeffH2Oice_Grundy1998] = load_abscoeffH2Oice_Grundy1998();


        wv = [abscoeffH2Oice_Grundy1998(26).data.wavelength]; % micro meter
        wv = wv(:)*1000; % convert to nano meter

        switch opt
            case 6
                dirpath_r = '/Users/itohy1/src/python/mie_yuki/';
                fnameList = {'Hapke_h2oIce_Grundy1998_150Kto220K.mat'};
                Hapke_H2Oicelib_Grundy1998 = [];
                for jdx = 1:length(fnameList)
                    fname = fnameList{jdx};
                    r = load(joinPath(dirpath_r,fname));
                    n = length(r.mD);
                    r_lib = struct('wavelength',cell(n,1),'r',cell(n,1),'libname',cell(n,1),...
                                   'unit_wavelength',cell(n,1),'name',cell(n,1));
                    for i=1:n
                        r_lib(i).wavelength = wv;
                        r_lib(i).r = r.r(:,i);
                        r_lib(i).unit_wavelength = 'nm';
                        r_lib(i).name = sprintf( ...
                            'R_H2Oice_Grundy_temp%3dK_mD%dum', ...
                            r.temp(i), r.mD(i) ...
                            );
                        r_lib(i).libname = fname;

                    end
                    Hapke_H2Oicelib_Grundy1998 = merge_struct(Hapke_H2Oicelib_Grundy1998,r_lib);
                end
            case 7
                dirpath_r = '/Users/itohy1/src/python/mie_yuki/';
                fnameList = {'Hapke_h2oIce_Grundy1998_150Kto270K.mat'};
                Hapke_H2Oicelib_Grundy1998 = [];
                for jdx = 1:length(fnameList)
                    fname = fnameList{jdx};
                    r = load(joinPath(dirpath_r,fname));
                    n = length(r.mD);
                    r_lib = struct('wavelength',cell(n,1),'r',cell(n,1),'libname',cell(n,1),...
                                   'unit_wavelength',cell(n,1),'name',cell(n,1));
                    for i=1:n
                        r_lib(i).wavelength = wv;
                        r_lib(i).r = r.r(:,i);
                        r_lib(i).unit_wavelength = 'nm';
                        r_lib(i).name = sprintf( ...
                            'R_H2Oice_Grundy_temp%3dK_mD%dum', ...
                            r.temp(i), r.mD(i) ...
                            );
                        r_lib(i).libname = fname;

                    end
                    Hapke_H2Oicelib_Grundy1998 = merge_struct(Hapke_H2Oicelib_Grundy1998,r_lib);
                end
            case 8
                dirpath_r = '/Users/itohy1/src/python/mie_yuki/';
                fnameList = {'Hapke_h2oIce_Grundy1998_150Kto270K_2.mat'};
                Hapke_H2Oicelib_Grundy1998 = [];
                for jdx = 1:length(fnameList)
                    fname = fnameList{jdx};
                    r = load(joinPath(dirpath_r,fname));
                    n = length(r.mD);
                    r_lib = struct('wavelength',cell(n,1),'r',cell(n,1),'libname',cell(n,1),...
                                   'unit_wavelength',cell(n,1),'name',cell(n,1));
                    for i=1:n
                        r_lib(i).wavelength = wv;
                        r_lib(i).r = r.r(:,i);
                        r_lib(i).unit_wavelength = 'nm';
                        r_lib(i).name = sprintf( ...
                            'R_H2Oice_Grundy_temp%3dK_mD%dum', ...
                            r.temp(i), r.mD(i) ...
                            );
                        r_lib(i).libname = fname;

                    end
                    Hapke_H2Oicelib_Grundy1998 = merge_struct(Hapke_H2Oicelib_Grundy1998,r_lib);
                end
            case 9
                dirpath_r = '/Users/itohy1/src/python/mie_yuki/';
                fnameList = {'Hapke_h2oIce_Grundy1998_150Kto270K_3.mat'};
                Hapke_H2Oicelib_Grundy1998 = [];
                for jdx = 1:length(fnameList)
                    fname = fnameList{jdx};
                    r = load(joinPath(dirpath_r,fname));
                    n = length(r.mD);
                    r_lib = struct('wavelength',cell(n,1),'r',cell(n,1),'libname',cell(n,1),...
                                   'unit_wavelength',cell(n,1),'name',cell(n,1));
                    for i=1:n
                        r_lib(i).wavelength = wv;
                        r_lib(i).r = r.r(:,i);
                        r_lib(i).unit_wavelength = 'nm';
                        r_lib(i).name = sprintf( ...
                            'R_H2Oice_Grundy_temp%3dK_mD%dum', ...
                            r.temp(i), r.mD(i) ...
                            );
                        r_lib(i).libname = fname;

                    end
                    Hapke_H2Oicelib_Grundy1998 = merge_struct(Hapke_H2Oicelib_Grundy1998,r_lib);
                end
            case 1
                dirpath_r = '/Users/itohy1/src/python/mie_yuki/';
                fnameList = {'Hapke_h2oIce_Grundy1998_150Kto230K_3.mat'};
                Hapke_H2Oicelib_Grundy1998 = [];
                for jdx = 1:length(fnameList)
                    fname = fnameList{jdx};
                    r = load(joinPath(dirpath_r,fname));
                    n = length(r.mD);
                    r_lib = struct('wavelength',cell(n,1),'r',cell(n,1),'libname',cell(n,1),...
                                   'unit_wavelength',cell(n,1),'name',cell(n,1));
                    for i=1:n
                        r_lib(i).wavelength = wv;
                        r_lib(i).r = r.r(:,i);
                        r_lib(i).unit_wavelength = 'nm';
                        r_lib(i).name = sprintf( ...
                            'R_H2Oice_Grundy_temp%3dK_mD%dum', ...
                            r.temp(i), r.mD(i) ...
                            );
                        r_lib(i).libname = fname;

                    end
                    Hapke_H2Oicelib_Grundy1998 = merge_struct(Hapke_H2Oicelib_Grundy1998,r_lib);
                end
            otherwise
                    error('opt %d is not defined',opt);
        end
        
        infoHapke_H2Oicelib_Grundy1998 = Hapke_H2Oicelib_Grundy1998;
        save(infocachefilepath,'infoHapke_H2Oicelib_Grundy1998');
    end
end

end