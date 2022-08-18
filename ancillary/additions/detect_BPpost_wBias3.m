function [BPpost_plus,band_bias_std,dev_sub,val_ratio_3d2,val_ratio_3d2_smth2] = detect_BPpost_wBias3(img,DMmask,varargin)
% [BPpost_plus,band_bias_std] = detect_BPpost_wBias2(BPall,img,DMmask,varargin)
% INPUTS
%  BPall: boolean, [1 x S x B]
%   all bad pixels
%  img: [L x S x B]
%   reference image
%  DMmask: boolean [1 x S x B]
% OU

bands = 1:252;
is_debug = false;
th_ratio = 0.005;
th_ratio_sigma = inf;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BANDS'
                bands = varargin{i+1};
            case 'DEBUG'
                is_debug = varargin{i+1};
            case 'THRESHOLD'
                th_ratio = varargin{i+1};
            case 'THRESHOLD_SIGMA'
                th_ratio_sigma = varargin{i+1};
            otherwise
                error('Unrecognized option %s',varargin{i});
        end
    end
end

if is_debug
    keyboard
end

% BPall = BPall==1;
% BPall = BPall(:,:,bands);
% 
% GPall = (BPall==0);
% GPall1nan = convertBoolTo1nan(GPall);

DMmask = DMmask==1;
DMmask = DMmask(:,:,bands);

% [L,S,B] = size(img);
% 
img_nanmed = nanmedian(img(:,:,bands),1);
% img_nanmed_GP1nan =  img_nanmed .* GPall1nan;
% img_nanmed_GP1nan(img_nanmed_GP1nan<0) = nan;

%%
% val_ratio_3d2 = img_nanmed_GP1nan./permute(img_nanmed_GP1nan,[3,2,1]);
val_ratio_3d2 = img_nanmed./permute(img_nanmed,[3,2,1]);

n=10;
AH = zeros(size(val_ratio_3d2,2),n);
x = reshape(1:size(val_ratio_3d2,2),[],1);
% x = x-median(x);
% x = x./max(abs(x))*n/2;
for j=1:n
    AH(:,j) = hermiteH(j-1,x);
end
AH = normalizevec(AH,1);

lambda_r = double(~isnan(val_ratio_3d2));

val_ratio_3d2_smth2 = nan(size(val_ratio_3d2));
for i=1:size(val_ratio_3d2_smth2,3)
    tic;
    clmns = DMmask(1,:,i);
    val_ratio_3d2_i = val_ratio_3d2(:,clmns,i);
    val_ratio_3d2_i(isnan(val_ratio_3d2_i)) = 1;
    
    [x,r,d,rho,Rhov,res_pv,res_dv,cost_val]...
        = lad_admm_gat_b(AH(clmns,:),val_ratio_3d2_i',...
                 'lambda_r',lambda_r(:,clmns,i)',...
                 'tol',1e-4,'maxiter',1000,'verbose','no',...
                 'PRECISION','double','gpu',false,'debug',false);
    % pmdl = AH*x;
    val_ratio_3d2_smth2(:,clmns,i) = (AH(clmns,:)*x)';
    toc;
end

r2 = val_ratio_3d2-val_ratio_3d2_smth2;
% r2 = val_ratio_3d2./val_ratio_3d2_smth2;
for i=1:size(r2,3)
    val_ratio_3d2(i,:,i) = nan;
    val_ratio_3d2_smth2(i,:,i) = nan;
    r2(i,:,i) = nan;
end
[dv_std] = robust_v3('med_abs_dev_from_med',r2,2,'data_center',0,'NOutliers',20);

% do computation again with regularization parameters.
val_ratio_3d2_smth2 = nan(size(val_ratio_3d2));
for i=1:size(val_ratio_3d2_smth2,3)
    tic;
    clmns = DMmask(1,:,i);
    val_ratio_3d2_i = val_ratio_3d2(:,clmns,i);
    val_ratio_3d2_i(isnan(val_ratio_3d2_i)) = 1;
    lambda_a = dv_std(:,:,i)'.*ones(size(AH,2),1);
    lambda_a([1,2],:) = 0;
    [x,r,d,rho,Rhov,res_pv,res_dv,cost_val]...
        = lad_admm_gat_b(AH(clmns,:),val_ratio_3d2_i',...
                 'lambda_r',lambda_r(:,clmns,i)',...
                 'lambda_a',lambda_a,...
                 'tol',1e-4,'maxiter',1000,'verbose','no',...
                 'PRECISION','double','gpu',false,'debug',false);
    % pmdl = AH*x;
    val_ratio_3d2_smth2(:,clmns,i) = (AH(clmns,:)*x)';
    toc;
end

r2 = val_ratio_3d2-val_ratio_3d2_smth2;
% r2 = val_ratio_3d2./val_ratio_3d2_smth2;
for i=1:size(r2,3)
    val_ratio_3d2(i,:,i) = nan;
    val_ratio_3d2_smth2(i,:,i) = nan;
    r2(i,:,i) = nan;
end

% Evaluate standard deviation of the profile in the cross-track direction.
% The standard deviations exhibit how each band is generally reliable.
[dv_std] = robust_v3('std',r2,2,'NOutliers',20);
dv_std = dv_std ./ nanmedian(val_ratio_3d2_smth2,2);
band_bias_r = nanmedian(dv_std,1);
% band_reliability(band_reliability>0.015) = inf;


% Next, evaluate how much each detector element is deviated. This is maded
% by taking the median 
[dv_ideal] = nanmedian(r2,1);

% evaluate likelihood
p = normpdf(dv_ideal,0,band_bias_r);
% p = normpdf(dv_ideal,1,band_reliability);

% bands proximity (closer bands are given larger weight)
b=4+squeeze(band_bias_r)./nanmedian(squeeze(band_bias_r));
band_proximity = 1./(2.*b).*exp(-1./b.*abs(permute(bands,[1,3,2]) - permute(bands,[2,1,3])));

% Weighted sum is used for estimating an ideal median values
r2_1nan = ones(size(r2));
r2_1nan(isnan(r2)) = nan;

img_nanmed_estimate = nansum(...
    val_ratio_3d2_smth2 .* permute(img_nanmed,[3,2,1]) ...
    .* band_proximity .* permute(p,[3,2,1])...
    ,1);

weight_sum = nansum(r2_1nan .* band_proximity  .* permute(p,[3,2,1]),1);

img_nanmed_estimate = img_nanmed_estimate ./ weight_sum;

% Evaluate the deviation from the ideal profile.
dev_sub = img_nanmed - img_nanmed_estimate;

% Expected band bias std
[band_bias_std] = robust_v3('std',dev_sub,2,'NOutliers',20);

% If the deviation is higher than th_ratio (default 0.005) or
% 10 sigma (default), the band is considered to be a bad pixel.
BPpost_plus = abs(dev_sub) > min(th_ratio,band_bias_std*th_ratio_sigma);


end