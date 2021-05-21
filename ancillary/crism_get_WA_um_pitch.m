function [WA_um_pitch] = get_WA_um_pitch_CRISM(WAdata)
% [WA_um_pitch] = get_WA_um_pitch_CRISM(WAdata)
%  Compute pitch of the wavelength frame
%  INPUTS
%  WAdata: CRISMdata obj, input wavelength frame
%  OUTPUTS
%  WA_um_pitch: pitch in micron of each channels defined by the data.

WA_nm = WAdata.readimgi(); % wavelength in nanometers
WA_nm(1,:,437) = nan;

WA_nm_extend = zeros([size(WA_nm,1),size(WA_nm,2),size(WA_nm,3)+2]);
WA_nm_extend(:,:,2:end-1) = WA_nm;
WA_nm_extend(:,:,1) = 2*WA_nm(:,:,1)-WA_nm(:,:,2);
WA_nm_extend(:,:,end)=2*WA_nm(:,:,end)-WA_nm(:,:,end-1);
WA_nm_between = (WA_nm_extend(:,:,2:end) + WA_nm_extend(:,:,1:end-1))/2;
WA_nm_pitch = WA_nm_between(:,:,2:end) - WA_nm_between(:,:,1:end-1);
WA_um_pitch = WA_nm_pitch .* (10^-3); % [um]

end
