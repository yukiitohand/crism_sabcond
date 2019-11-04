function [out] = sabcondv3_pub_column_test(obs_id,c,varargin)
% [out] = sabcondv3_pub_column_test(obs_id,columns,varargin)
%   test sabcondv3_pub only at one column.
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
%   see sabcondv3_pub
%   **Default values different from sabcondv3_pub**
%       'INTERLEAVE_OUT'    : 'bls'
%       'SUBSET_COLUMNS_OUT': true
%       'SAVE_FILE'         : false

[out] = sabcondv3_pub(obs_id,'INTERLEAVE_OUT','bls',...
    'SUBSET_COLUMNS_OUT',1,'SAVE_FILE',false, varargin{:},'Columns',c);

end