function [] = crism_sabcond_startup_addpath()
%-------------------------------------------------------------------------%
% % Automatically find the path to toolboxes
fpath_self = mfilename('fullpath');
[dirpath_self,filename] = fileparts(fpath_self);
mtch = regexpi(dirpath_self,'(?<parent_dirpath>.*)/crism_sabcond[/]{0,1}','names');
toolbox_root_dir = mtch.parent_dirpath;

%-------------------------------------------------------------------------%
% name of the directory of each toolbox
base_toolbox_dirname  = 'base';
envi_toolbox_dirname  = 'envi';
pds3_toolbox_dirname  = 'pds3_toolbox';
crism_toolbox_dirname = 'crism_toolbox';
crmcal_toolbox_dirname        = 'crism_calibration';
crism_sabcond_toolbox_dirname = 'crism_sabcond';
cntrmvl_toolbox_dirname       = 'continuumRemoval';
%-------------------------------------------------------------------------%
pathCell = strsplit(path, pathsep);

%% base toolbox
base_toolbox_dir = [toolbox_root_dir '/' base_toolbox_dirname ];
% joinPath in base toolbox will be used in the following. "base" toolbox
% need to be loaded first. base/joinPath.m automatically determine the
% presence of trailing slash, so you do not need to worry it.
if exist(base_toolbox_dir,'dir')
    if ~check_path_exist(base_toolbox_dir, pathCell)
        addpath(base_toolbox_dir);
    end
else
    warning([ ...
        'base toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/base/'
        ]);
end

%% envi toolbox
envi_toolbox_dir = joinPath(toolbox_root_dir, envi_toolbox_dirname);

if exist(envi_toolbox_dir,'dir')
    if ~check_path_exist(envi_toolbox_dir, pathCell)
        run(joinPath(envi_toolbox_dir,'envi_startup_path'));
    end
else
    warning([ ...
        'envi toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/envi/'
        ]);
end

%% pds3_toolbox
pds3_toolbox_dir = joinPath(toolbox_root_dir, pds3_toolbox_dirname);

if exist(pds3_toolbox_dir,'dir')
    if ~check_path_exist(pds3_toolbox_dir, pathCell)
        run(joinPath(pds3_toolbox_dir,'pds3_startup_path'));
    end
else
    warning([ ...
        'pds3_toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/pds3_toolbox/'
        ]);
end

%% crism_toolbox
crism_toolbox_dir = joinPath(toolbox_root_dir, crism_toolbox_dirname);

if exist(crism_toolbox_dir,'dir')
    if ~check_path_exist(crism_toolbox_dir, pathCell)
        run(joinPath(crism_toolbox_dir,'crism_startup_path'));
    end
else
    warning([ ...
        'crism_toolbox is not detected. Download from' '\n' ...
        '   https://github.com/yukiitohand/crism_toolbox/'
        ]);
end

%%
% future implementation
% continuum_removal toolbox
% crism_calibration

%% crism_sabcond toolbox
crism_sabcond_toolbox_dir = joinPath(toolbox_root_dir, crism_sabcond_toolbox_dirname);
if ~check_path_exist(pds3_toolbox_dir, pathCell)
    addpath( ...
        crism_sabcond_toolbox_dir                                  , ...
        joinPath(crism_sabcond_toolbox_dir,'legacy/')              , ...
        joinPath(crism_sabcond_toolbox_dir,'legacy/sabcondc/')     , ...
        joinPath(crism_sabcond_toolbox_dir,'ancillary/')           , ...
        joinPath(crism_sabcond_toolbox_dir,'ancillary/additions/') , ...
        joinPath(crism_sabcond_toolbox_dir,'ancillary/library/')   , ...
        joinPath(crism_sabcond_toolbox_dir,'ancillary/T/')         , ...
        joinPath(crism_sabcond_toolbox_dir,'ancillary/util/')      , ...
        joinPath(crism_sabcond_toolbox_dir,'sabcondc/')              ...
    );
end

end

%%
function [onPath] = check_path_exist(dirpath, pathCell)
    % pathCell = strsplit(path, pathsep, 'split');
    if ispc || ismac 
      onPath = any(strcmpi(dirpath, pathCell));
    else
      onPath = any(strcmp(dirpath, pathCell));
    end
end