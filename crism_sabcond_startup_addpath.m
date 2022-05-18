function [] = crism_sabcond_startup_addpath()
%% Automatically find the path to toolboxes
fpath_self = mfilename('fullpath');
[dirpath_self,filename] = fileparts(fpath_self);
path_ptrn = '(?<parent_dirpath>.*)/(?<toolbox_dirname>[^/]+[/]{0,1})';
mtch = regexpi(dirpath_self,path_ptrn,'names');
toolbox_root_dir = mtch.parent_dirpath;
crism_sabcond_toolbox_dirname  = mtch.toolbox_dirname;

%% find dependent toolboxes
dList = dir(toolbox_root_dir);
error_if_not_unique = true;
pathCell = strsplit(path, pathsep);

% base toolbox
[base_toolbox_dir,base_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'base',error_if_not_unique);
if ~check_path_exist(base_toolbox_dir, pathCell)
    addpath(base_toolbox_dir);
end

% envi toolbox
[envi_toolbox_dir,envi_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'envi',error_if_not_unique);
if ~check_path_exist(envi_toolbox_dir, pathCell)
    run(joinPath(envi_toolbox_dir,'envi_startup_addpath'));
end

% pds3_toolbox
[pds3_toolbox_dir,pds3_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'pds3_toolbox',error_if_not_unique);
if ~check_path_exist(pds3_toolbox_dir, pathCell)
    run(joinPath(pds3_toolbox_dir,'pds3_startup_addpath'));
end

% suwab toolbox
[suwab_toolbox_dir,suwab_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'suwab_toolbox',error_if_not_unique);
if ~check_path_exist(suwab_toolbox_dir, pathCell)
    run(joinPath(suwab_toolbox_dir,'suwab_toolbox_startup_addpath'));
end

% continuumRemoval toolbox
[cntrmvl_toolbox_dir,cntrmvl_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'continuumRemoval',error_if_not_unique);
if ~check_path_exist(cntrmvl_toolbox_dir, pathCell)
    run(joinPath(cntrmvl_toolbox_dir,'cntrmvl_toolbox_startup_addpath'));
end

% crism_toolbox
[crism_toolbox_dir,crism_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'crism_toolbox',error_if_not_unique);
if ~check_path_exist(crism_toolbox_dir, pathCell)
    run(joinPath(crism_toolbox_dir,'crism_startup_addpath'));
end

% crism_calibration toolbox
[crmcal_toolbox_dir,crmcal_toolbox_dirname] = get_toolbox_dirname(dList, ...
    'crism_calibration',error_if_not_unique);
if ~check_path_exist(crmcal_toolbox_dir, pathCell)
    run(joinPath(crmcal_toolbox_dir,'crism_calibration_startup_addpath'));
end

%% crism_sabcond toolbox
crism_sabcond_toolbox_dir = joinPath(toolbox_root_dir, crism_sabcond_toolbox_dirname);
if ~check_path_exist(crism_sabcond_toolbox_dir, pathCell)
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

function [toolbox_dirpath,toolbox_dirname] = get_toolbox_dirname(dList,toolbox_dirname_wover,error_if_not_unique)
    dirname_ptrn = sprintf('(?<toolbox_dirname>%s(-[\\d\\.]+){0,1}[/]{0,1})',toolbox_dirname_wover);
    mtch_toolbox_dirname = regexpi({dList.name},dirname_ptrn,'names');
    mtchidx = find(not(cellfun('isempty',mtch_toolbox_dirname)));
    toolbox_root_dir = dList(1).folder;
    if length(mtchidx)==1
        toolbox_dirname = mtch_toolbox_dirname{mtchidx}.toolbox_dirname;
        toolbox_dirpath = [toolbox_root_dir, '/', toolbox_dirname];
    elseif isempty(mtchidx)
        toolbox_dirname = [];
        toolbox_dirpath = [];
        if error_if_not_unique
            
            error(['Cannot find %s toolbox in\n %s', ...
                   'Download from\n', ...
                   '   https://github.com/yukiitohand/%s/', ...
                  ], ...
                toolbox_dirname_wover,toolbox_root_dir,toolbox_dirname_wover);
        end
    else % length(mtchidx)>1
        toolbox_dirname = {cat(2,mtch_toolbox_dirname{mtchidx}).toolbox_dirname};
        toolbox_dirpath = cellfun(@(x) strjoin({toolbox_root_dir,x},'/'), toolbox_dirname, ...
            'UniformOutput',false);
        if error_if_not_unique
            toolbox_root_dir = dList(1).folder;
            error('Multiple %s toolboxes %s are detected in\n %s', ...
            toolbox_dirname_wover,strjoin(toolbox_dirname,','), toolbox_root_dir);
        end
        
    end
end