function [Yim_1ord] = approx_HSI3d_1ord(Yim,WA,varargin)
% [Yim_1ord] = approx_HSI3d_1ord(Y,WA,varargin)
%  approximate each column of hyperspectral image with a first order matrix. 
%  First order matrix can be estimated either by
%    1) using a least square on the column mean. 
%    2) using SVD
%  Yim: [L x S x B] hyeprspectral image
%  x: [1 x S x B], wavelength frame
%  
%  OUTPUTS
%  Yim_1ord: [L x S x B] approximated image cube
%  OPTIONAL PARAMETERS
%  'MODE' : mode for approximation {'COLUMN_MEAN','COLUMN_MEDIAN','COLUMN_MED_L1','SVD'}
%           (default) 'COLUMN_MEAN'
%  'BAD': [L x S x B] Boolean 3d array, 1 if the value should be replaced 0 if
%         not. In any case, NaN
%         will be replaced.
%         (default) false([L,S,B])

[L,S,B] = size(Yim);
Yim_bad = false([L,S,B]);
app_mode = 'COLUMN_MEAN';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BAD'
                Yim_bad = varargin{i+1};
                validateattributes(Yim_bad,'logical',{'size',[L,S,B]});
            case 'MODE'
                app_mode = varargin{i+1};
                app_mode = validatestring(upper(app_mode),...
                    {'COLUMN_MEAN','COLUMN_MEDIAN','COLUMN_MED_L1','SVD'},...
                    mfilename,'MODE');
            otherwise
                % Hmmm, something wrong with the parameter string
                error('Unrecognized option: %s',varargin{i});   
        end
    end
end

Yim = permute(Yim,[3,1,2]);
WA = permute(WA,[3,1,2]);
Yim_1ord = nan([B,L,S]);

Yim_bad = permute(Yim_bad,[3,1,2]);
nisnanWA = ~isnan(WA);

for s=1:S
    if any(nisnanWA(:,:,s))
        [Yim_1ord(nisnanWA(:,:,s),:,s)] = approx_Mat_1ord(...
            Yim(nisnanWA(:,:,s),:,s),WA(nisnanWA(:,:,s),:,s),...
            'MODE',app_mode,'BAD',Yim_bad(nisnanWA(:,:,s),:,s));
    end
end

Yim_1ord = permute(Yim_1ord,[2,3,1]);

end