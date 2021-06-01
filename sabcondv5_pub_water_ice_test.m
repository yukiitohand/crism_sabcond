function [result] = sabcondv5_pub_water_ice_test(obs_id,opt_icelib,varargin)
% [out] = sabcondv5_pub_water_ice_test(obs_id,varargin)
%   test the presence of water ice using sabcondv3_pub.
%
% INPUT Parameters
%   obs_id : string, 
%       observation ID of the image to be processed. Currently, EPF 
%       measurements accompanied with scene measurements are not supported.
%   opt_icelib: integer,
%       index id for ice library to be used.
% OUTPUT Parameters
%   result: struct, storing result of the test.
%     *fields* (C indicates column, L indicates lines, )
%      'presence_H2Oice' : boolean, 
%           true if test is true and false otherwise
%      'presence_H2Oice_columns' : boolean array [1 C]
%           whether or not Xice_column_mean_bland is greater than 
%           THRESHOLD_ICE_ABUNDANCE
%      'Xice_column_mean_bland' : array [1,C], 
%           column mean of H2O ice abundances of bland spectra
%      'Xice_map' : 2d array [L C], 
%           map of total H2O ice abundances
%      'map_bland' : 2d boolean array [L C], 
%           map of bland spectra
%      'non_bland_mode: boolean,
%           whether or not only bland spectra are used for test
%           true if all spectra are used, false only bland spectra are used
%
%      'columns' : integer array, column indexes used for this test
%      'bands'   : integer array, bands used for processing
%      'lines'   : integer array, lines used for processing
%
% OPTIONAL Parameters
%   'THRESHOLD_BLAND': scalar
%       threshold value for which spectra are considered bland.
%       (default) 0.5
%   'THRESHOLD_ICE_ABUNDANCE': sclar
%       threshold value for abundance of ice.
%       (default) 0.1
%   'THRESHOLD_ICE_DOMINANCE': sclar
%       threshold value for spatial dominance of ice
%       (default) 0.8
%     
%   see sabcondv5_pub for rest of the parameters.
%   **Default values different from sabcondv3_pub**
%       'SUBSET_COLUMNS_OUT': true
%       'SAVE_FILE'         : false
%       'ALIB_OUT'          : true
%       'nIter'             : 0
%       'COLUMN_SKIP'       : 50

threshold_bland = 0.5;
threshold_Xice  = 0.1;
threshold_ice_dominance = 0.8;
column_skip = 50;

idx_varargin_rmvl = [];

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            % ## THRESHOLDS #----------------------------------------------
            case 'THREHOLD_BLAND'
                threshold_bland = varargin{i+1};
                idx_varargin_rmvl = [idx_varargin_rmvl i i+1];
            case 'THRESHOLD_ICE_ABUNDANCE'
                threshold_Xice = varargin{i+1};
                idx_varargin_rmvl = [idx_varargin_rmvl i i+1];
            case 'THREHOLD_ICE_DOMINANCE'
                threshold_ice_dominance = varargin{i+1};
                idx_varargin_rmvl = [idx_varargin_rmvl i i+1];
        end
    end
end

idx_varargin_retained = setdiff(1:length(varargin),idx_varargin_rmvl);
varargin = varargin(idx_varargin_retained);

[out] = sabcondv5_pub(obs_id,'OPT_ICELIB',opt_icelib,...
    'INTERLEAVE_OUT','lsb','SUBSET_COLUMNS_OUT',1,'SAVE_FILE',false,...
    'Alib_out', true,'nIter',0,'COLUMN_SKIP',column_skip,varargin{:});

% select bland spectra
map_bland = abs(sum(log(out.AB_est(:,:,out.bands)),3)) < threshold_bland;
%
map_bland_1nan = convertBoolTo1nan(map_bland);
% amount of ice
Xice_map = nan(length(out.lines),length(out.columns));
for i=1:length(out.columns)
    if ~isempty(out.ancillaries(i).Xice)
        Xice_map(:,i) = sum(out.ancillaries(i).Xice,1);
    end
end
Xice_column_mean_bland = mean(Xice_map.*map_bland_1nan,1,'omitnan');
columns_not_consider = all(isnan(map_bland_1nan),1);
if all(columns_not_consider)
    fprintf('No bland spectra are detected on selected columns.\n');
    fprintf('Test will be performed on non-bland spectra\n');
    non_bland_mode = true;
    Xice_column_mean_bland = mean(Xice_map,1,'omitnan');
    columns_not_consider = all(isnan(Xice_map),1);
else
    non_bland_mode = false;
end

Xice_column_mean_bland(columns_not_consider) = nan;

presence_H2Oice_columns = double(Xice_column_mean_bland > threshold_Xice);
presence_H2Oice_columns(columns_not_consider) = nan;

presence_H2Oice = nanmean(presence_H2Oice_columns,2) > threshold_ice_dominance;

result = [];
result.presence_H2Oice         = presence_H2Oice;
result.presence_H2Oice_columns = presence_H2Oice_columns;
result.Xice_column_mean_bland  = Xice_column_mean_bland;
result.Xice_map         = Xice_map;
result.map_bland        = map_bland;
result.columns          = out.columns;
result.bands            = out.bands;
result.lines            = out.lines;
result.non_bland_mode   = non_bland_mode;


end