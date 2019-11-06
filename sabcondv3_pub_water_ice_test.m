function [result] = sabcondv3_pub_water_ice_test(obs_id,opt_icelib,varargin)
% [out] = sabcondv3_pub_water_ice_test(obs_id,varargin)
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
%     *fields*
%      water_ice_result:
%      water_ice_exist:
%      Xice_mean_insig:
%      Xice_map:
%      insig_map:
%      columns:
%      bands:
%      lines:
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
%   see sabcondv3_pub for rest of the parameters.
%   **Default values different from sabcondv3_pub**
%       'SUBSET_COLUMNS_OUT': true
%       'SAVE_FILE'         : false
%       'ALIB_OUT'          : true
%       'nIter'             : 0
%       'COLUMN_SKIP'       : 50

threshold_insig_logAB_bprmvd = 0.5;
threshold_Xice = 0.1;
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
                threshold_insig_logAB_bprmvd = varargin{i+1};
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

idx_varargin_retained = setdiff(length(varargin),idx_varargin_rmvl);
varargin = varargin(idx_varargin_retained);

[out] = sabcondv3_pub(obs_id,'OPT_ICELIB',opt_icelib,...
    'INTERLEAVE_OUT','lsb','SUBSET_COLUMNS_OUT',1,'SAVE_FILE',false,...
    'Alib_out', true,'nIter',0,'COLUMN_SKIP',column_skip,varargin{:});

% select bland spectra
idx_insig = abs(sum(log(out.AB_est(:,:,out.bands),3))) ...
                                     < threshold_insig_logAB_bprmvd;

% amount of ice
Xice = cat(3,out.ancillaries.Xice);
Xice = reshape(Xice,[2,3,1]);
Xice_map = sum(Xice,3);
Xice_mean_insig = mean(Xice_map(idx_insig));

water_ice_exist = Xice_mean_insig > threshold_Xice;

if nanmean(water_ice_exist) > threshold_ice_dominance
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