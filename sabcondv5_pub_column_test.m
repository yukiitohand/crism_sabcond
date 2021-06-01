function [out] = sabcondv5_pub_column_test(obs_id,c,varargin)
% [out] = sabcondv5_pub_column_test(obs_id,columns,varargin)
%   test sabcondv5_pub only at one column.
%
% INPUT Parameters
%   obs_id : string, 
%       observation ID of the image to be processed. Currently, EPF 
%       measurements accompanied with scene measurements are not supported.
%   c: array or scalar
%       column to be processed
% OUTPUT Parameters
%   see sabcondv3_pub
%
% OPTIONAL Parameters
%   see sabcondv5_pub
%   **Default values different from sabcondv3_pub**
%       'INTERLEAVE_OUT'    : 'bls'
%       'SUBSET_COLUMNS_OUT': true
%       'SAVE_FILE'         : false
%       'ALIB_OUT'          : true

[out] = sabcondv5_pub(obs_id,'INTERLEAVE_OUT','bls',...
    'SUBSET_COLUMNS_OUT',1,'SAVE_FILE',false, 'Alib_out', true,...
    varargin{:},'Columns',c);

end