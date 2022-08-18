function [Y_1ord] = approx_Mat_1ord(Y,x,varargin)
% [Y_1ord] = approx_Mat_1ord(Y,x,varargin)
%  approximate Y with a first order matrix. First order matrix can be 
%  estimated either by
%    1) using a least square on the column mean. 
%    2) using SVD
%  Y: a matrix, [M x N]
%  x: [M x 1] vector, sample points for each row of Y.
%  
%  OUTPUTS
%  Y_1ord: [M x N] approximated matrix
%  OPTIONAL PARAMETERS
%  'MODE' : mode for approximation {'COLUMN_MEAN','COLUMN_MEDIAN','COLUMN_MED_L1','SVD'}
%           (default) 'COLUMN_MEAN'
%  'BAD': [M x N] Boolean matrix, 1 if the value should be replaced 0 if
%         not. If empty, nothing is considered as bad. In any case, NaN
%         will be replaced.
%         (default) []

Y_bad = [];
app_mode = 'COLUMN_MEAN';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BAD'
                Y_bad = varargin{i+1};
            case 'MODE'
                app_mode = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error('Unrecognized option: %s',varargin{i});   
        end
    end
end

if any(isnan(x))
    error('no elements of x can be NaN.');
end

% replace invalid values (NaN) with interpolation
Y_isnan = isnan(Y);
if ~isempty(Y_bad)
    Y_bad = or(Y_bad,Y_isnan);
else
    Y_bad = Y_isnan;
end

Y_filled = interp_nan_column(Y,Y_bad,x);

switch upper(app_mode)
    case 'COLUMN_MEAN'
        % column mean
        mYifcol = nanmean(Y_filled,2);
        % row component estimation by least square
        mYifrow = mYifcol'*Y_filled/norm(mYifcol,2)^2;
        % approximation
        Y_1ord = mYifcol*mYifrow;
    case 'COLUMN_MED'
        % column mean
        mYifcol = nanmedian(Y_filled,2);
        % row component estimation by least square
        mYifrow = mYifcol'*Y_filled/norm(mYifcol,2)^2;
        % approximation
        Y_1ord = mYifcol*mYifrow;
    case 'SVD'
        [U,S,V] = svds(Y_filled,1);
        Y_1ord = U * S * V';
    case 'COLUMN_MED_L1'
        % column mean
        mYifcol = nanmedian(Y_filled,2);
        % row component estimation by least square
        mYifrow = wlad_gadmm_a_v2(mYifcol,Y_filled,'tol',1e-10,'maxiter',1000,'verbose',false);
        % approximation
        Y_1ord = mYifcol*mYifrow;
    otherwise
        error('Undefined mode');
end

end