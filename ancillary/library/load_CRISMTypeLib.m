function [crismTypeLib] = load_CRISMTypeLib(opt_CRISMTypeLib)
% load MRO CRISM Type Spectra Library
% Input Parameters
%   opt_CRISMTypeLib : subclass id
% Output Parameters
%   crismTypeLib: struct, see 'readCRISMTypeLibrary.m' for detail
%         additional field
%            'spclib': 'CRISM Type Spectral Library'
%            'cumidex'
%            'reflectannce'
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
    otherwise
        error('opt_CRISMTypeLib %d is not defined',opt_CRISMTypeLib);
end

if ~isempty(crismTypeLib)
    [crismTypeLib.spclib] = deal('CRISM Type Spectral Library');
    iList = num2cell(1:length(crismTypeLib));
    [crismTypeLib.cumindex] = iList{:};
end

end