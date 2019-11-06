function [result] = sabcondv3_pub_water_ice_test(obs_id,varargin)
% [out] = sabcondv3_pub_water_ice_test(obs_id,varargin)
%   test the presence of water ice using sabcondv3_pub.
%
% INPUT Parameters
%   obs_id : string, 
%       observation ID of the image to be processed. Currently, EPF 
%       measurements accompanied with scene measurements are not supported.
% OUTPUT Parameters
%   see sabcondv3_pub
%
% OPTIONAL Parameters
%   see sabcondv3_pub
%   **Default values different from sabcondv3_pub**
%       'SUBSET_COLUMNS_OUT': true
%       'SAVE_FILE'         : false
%       'ALIB_OUT'          : true
%       'nIter'             : 0

threshold_insig_logAB_bprmvd = 0.5;
threshold_Xice = 0.1;
column_skip = 50;

[out] = sabcondv3_pub(obs_id,'INTERLEAVE_OUT','lsb',...
    'SUBSET_COLUMNS_OUT',1,'SAVE_FILE',false, 'Alib_out', true,...
    'nIter',0,'COLUMN_SKIP',column_skip,varargin{:});

% select bland spectra
idx_insig = abs(sum(log(out.AB_est(:,:,out.bands),3))) ...
                                     < threshold_insig_logAB_bprmvd;

% amount of ice
Xice_map = sum(out.ancillariesXice,3);
Xice_mean_insig = mean(Xice_map(idx_insig));

water_ice_exist = Xice_mean_insig > threshold_Xice;

if nanmean(water_ice_exist)>0.8
    water_ice_result = 1;
else
    water_ice_result = 0;
end


result = [];
result.water_ice_result = water_ice_result;
result.water_ice_exist  = water_ice_exist;
result.Xice_mean_insig  = Xice_mean_insig;
result.Xice_map         = Xice_map;
result.insig_map        = idx_insig;
result.columns          = out.columns;
result.bands            = out.bands;
result.lines            = out.lines;


end