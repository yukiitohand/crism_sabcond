function [Yim_filled] = interp_nan_column_3D(Yim,WA,varargin)
% [Yim_filled] = interp_nan_column_3D(Yim,WA,varargin)
%  3d version of interp_nan_column
%  INPUTS
%  Yim: [L x S x B] hyeprspectral image
%  WA: [1 x S x B], wavelength frame
%  
%  OUTPUTS
%  Yim_1ord: [L x S x B] approximated image cube
%  OPTIONAL PARAMETERS
%  'BAD': [L x S x B] Boolean 3d array, 1 if the value should be replaced 0 if
%         not. In any case, NaN
%         will be replaced.
%         (default) false([L,S,B])

[L,S,B] = size(Yim);
Yim_bad = false([L,S,B]);

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BAD'
                Yim_bad = varargin{i+1};
                validateattributes(Yim_bad,'logical',{'size',[L,S,B]});
            otherwise
                % Hmmm, something wrong with the parameter string
                error('Unrecognized option: %s',varargin{i});   
        end
    end
end

% replace invalid values (NaN) with interpolation
Yim_isnan = isnan(Yim);
if ~isempty(Yim_bad)
    Yim_bad = or(Yim_bad,Yim_isnan);
else
    Yim_bad = Yim_isnan;
end

Yim = permute(Yim,[3,1,2]);
WA = permute(WA,[3,1,2]);
Yim_filled = nan([B,L,S]);
Yim_bad = permute(Yim_bad,[3,1,2]);
nisnanWA = ~isnan(WA);

for s=1:S
    if any(nisnanWA(:,:,s))
        [Yim_filled(nisnanWA(:,:,s),:,s)] = interp_nan_column(...
            Yim(nisnanWA(:,:,s),:,s),Yim_bad(nisnanWA(:,:,s),:,s),...
            WA(nisnanWA(:,:,s),:,s));
    end
end

Yim_filled = permute(Yim_filled,[2,3,1]);

end