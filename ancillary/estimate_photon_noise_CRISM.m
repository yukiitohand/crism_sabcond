function [photon_noise_mad_stdif,WA_um_pitch] = estimate_photon_noise_CRISM(TRRYRAdata)
% [photon_noise_mad_stdif] = estimate_photon_noise_CRISM(TRRYRAdata)
%  Estimate the median absolute deviation of photon noise in the domain of
%  I/F. 
%  INPUTS
%  TRRYRAdara: input radiance image.
%  OUTPUTS
%  photon_noise_mad_stdif: the median absolute deviation of photon noise 
%     in the domain of I/F. [L x S x B]
%  WA_um_pitch: 

if isempty(TRRYRAdata.basenamesCDR)
    TRRYRAdata.load_basenamesCDR();
end

RDimg = TRRYRAdata.readimgi();
WAdata = TRRYRAdata.readCDR('WA');
WA_nm = WAdata.readimgi(); % wavelength in nanometers
WA_nm(1,:,437) = nan;

WA_nm_extend = zeros([size(WA_nm,1),size(WA_nm,2),size(WA_nm,3)+2]);
WA_nm_extend(:,:,2:end-1) = WA_nm;
WA_nm_extend(:,:,1) = 2*WA_nm(:,:,1)-WA_nm(:,:,2);
WA_nm_extend(:,:,end)=2*WA_nm(:,:,end)-WA_nm(:,:,end-1);
WA_nm_between = (WA_nm_extend(:,:,2:end) + WA_nm_extend(:,:,1:end-1))/2;
WA_nm_pitch = WA_nm_between(:,:,2:end) - WA_nm_between(:,:,1:end-1);
WA_um_pitch = WA_nm_pitch .* (10^-3); % [um]

SFdata = TRRYRAdata.readCDR('SF');
SFimg = SFdata.readimgi();

% [RDimg_1ord] = approx_HSI3d_1ord(RDimg,WA_nm,'MODE','COLUMN_MEAN');
[RDimg_filled] = interp_nan_column_3D(RDimg,WA_nm);

[photon_noise_mad_stdif] = estimate_photon_noise_CRISM_base(RDimg_filled,WA_nm,WA_um_pitch,TRRYRAdata.lbl,SFimg);

end
