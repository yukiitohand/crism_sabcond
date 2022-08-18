function [dev_coef,img_nanmed_estimate,dev_coef_ori] = calc_deviation_BPpost2(BP,img)
% [dev_coef] = calc_deviation_BPpost(BPpri,BPall,img)
% INPUTS
%  BPpri: boolean, [1 x S x B]
%   bad pixels detected at the early stage of the bad pixel detection
%  BPall" boolean, [1 x S x B]
%   all bad pixels
%  img: [L x S x B]
%   reference image

BP = BP==1;

[L,S,B] = size(img);

img_nanmed = median(img,1,'omitnan');

GP = (BP==0);
GP1nan = convertBoolTo1nan(GP);

img_nanmed_GP1nan =  img_nanmed .* GP1nan;
img_nanmed_GP1nan(img_nanmed_GP1nan<0) = nan;

n_old = sum(isnan(img_nanmed_GP1nan),'all');
n_new = 1;

while n_new<n_old
    n_old = sum(isnan(img_nanmed_GP1nan),'all');
    %fprintf('%d,%d,%d\n',sum(isnan(img_nanmed),'all'),sum(isnan(img_nanmed_GP1nan),'all'));
    img_nanmed_GP1nan_isnan = or(isnan(img_nanmed_GP1nan),img_nanmed_GP1nan<0);
    
    %-------------------------------------------------------------------------%
    % first take the ratios of the band neighbors (left and right separately)
    val_ratio_bl = nan(1,S,B);
    val_ratio_bl(:,:,2:B) = img_nanmed_GP1nan(:,:,2:B)./img_nanmed_GP1nan(:,:,1:(B-1));
    val_ratio_br = nan(1,S,B);
    val_ratio_br(:,:,1:(B-1)) = img_nanmed_GP1nan(:,:,1:(B-1))./img_nanmed_GP1nan(:,:,2:B);
    
    %-------------------------------------------------------------------------%
    % first take the ratios of the further band neighbors (left and right separately)
    % val_ratio_bl2 = nan(1,S,B);
    % val_ratio_bl2(:,:,3:B) = img_nanmed_GP1nan(:,:,3:B)./img_nanmed_GP1nan(:,:,1:(B-2));
    % val_ratio_br2 = nan(1,S,B);
    % val_ratio_br2(:,:,1:(B-2)) = img_nanmed_GP1nan(:,:,1:(B-2))./img_nanmed_GP1nan(:,:,3:B);

    %-------------------------------------------------------------------------%
    % second, averaging the ratios among the spatial neighbors
    val_ratio_sbl = nan(1,S,B);
    val_ratio_sbl(:,2:(S-1),:) = mean(cat(1,val_ratio_bl(:,1:(S-2),:),val_ratio_bl(:,3:S,:)),1,'omitnan');

    val_ratio_sbr = nan(1,S,B);
    val_ratio_sbr(:,2:(S-1),:) = mean(cat(1,val_ratio_br(:,1:(S-2),:),val_ratio_br(:,3:S,:)),1,'omitnan');

%     val_ratio_combined = nan(1,S,B);
%     a = cat(1,...
%         val_ratio_sbl(:,:,2:(B-1)).*img_nanmed_GP1nan(:,:,1:(B-2)),...
%         val_ratio_sbr(:,:,2:(B-1)).*img_nanmed_GP1nan(:,:,3:B));
%     a(a<0) = nan;
%     val_ratio_combined(:,:,2:(B-1)) = geomean(a,1,'omitnan');
    
    val_ratio_combined = nan(2,S,B);
    val_ratio_combined(1,:,2:B) = val_ratio_sbl(:,:,2:B).*img_nanmed_GP1nan(:,:,1:(B-1));
    val_ratio_combined(2,:,1:(B-1)) = val_ratio_sbr(:,:,1:(B-1)).*img_nanmed_GP1nan(:,:,2:B);
    val_ratio_combined(val_ratio_combined<0) = nan;
    val_ratio_combined = geomean(val_ratio_combined,1,'omitnan');
    
    img_nanmed_GP1nan(img_nanmed_GP1nan_isnan) = val_ratio_combined(img_nanmed_GP1nan_isnan);
    
    % aaa = cat(1,aaa,val_ratio_combined);
    
    n_new = sum(isnan(img_nanmed_GP1nan),'all');
end

dev_coef = img_nanmed./img_nanmed_GP1nan;

dev_coef_ori = dev_coef;

dev_coef(GP) = 1;

img_nanmed_estimate = img_nanmed_GP1nan;

end



