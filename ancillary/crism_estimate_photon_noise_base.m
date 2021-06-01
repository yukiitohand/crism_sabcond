function [photon_noise_mad_stdif] = crism_estimate_photon_noise_base(...
    RDimg,WA_nm,WA_um_pitch,lbl,SFimg)
% [photon_noise_mad_stdif] = crism_estimate_photon_noise_base(...
%    RDimg,WA_nm,WA_um_pitch,lbl,SFimg)
%  Estimate the median absolute deviation of photon noise in the domain of
%  I/F. 
%  INPUTS
%  RDimg: input radiance image.
%  WA_nm: wavelength frame
%  WA_um_pitch: width for each wavelength
%  lbl: CRISM LBL file storing acquisition information
%  SFimg: Irradiance image 
%  *Note:
%    RDimg, WA_nm, WA_um_pitch, and SFimg need to have sizes compatible 
%    with implicit expansion.
%    
%  OUTPUTS
%  photon_noise_mad_stdif: the median absolute deviation of photon noise 
%     in the domain of I/F.


%% Some physical parameters
h_planck = 6.626070040 * 10^(-34); % planck constant [kg*m^2/s]
c_light = 299792458; % [m/s] speed of light
ap       = 0.1; %[m], aperture 100mm
alt      = 3.0*10^5; %[m], altitude of the sattelite 300km
fl       = 0.441; %[m], focal length 441mm
slw      = 27*10^(-6); %[m] slit width 27um
slh_1del = 27*10^(-6);%[m] 1 detector element equivalent slit height 27um
eta = 0.8; % Optical and quantum efficiency (rough guess...)

%%

%1 pixel equivalent slit height, considering spatial binning
binx = lbl.PIXEL_AVERAGING_WIDTH;
slh_1pxl = slh_1del * binx; % [m]

%% Computation of the statistic of photon noise
sr = pi*((ap/2).^2) / (alt)^2; % [steradians] solid angle
ifov_along_track = alt/fl * slw;
ifov_cross_track = alt/fl * slh_1pxl;
A = ifov_cross_track*ifov_along_track; %[m^2] area from which the light is comming

% get integration time
integ_t = lbl.MRO_EXPOSURE_PARAMETER;
rateHz = lbl.MRO_FRAME_RATE.value;
[tms] = crism_get_integrationTime(integ_t,rateHz,'Hz');
tsec = tms/1000; %[s] integration time 


WA_m = WA_nm * (10.^(-9)); % [m] converting to meters

d_km = lbl.SOLAR_DISTANCE.value;
[ d_au ] = km2au( d_km );

% photon_cts = RDimg .* (A*sr*tsec).* WA_um_pitch ./ ((h_planck.*c_light)./WA_m) .* eta;
photon_cts = (eta.*A*sr*tsec./ (h_planck.*c_light) ) .*  (WA_um_pitch.*WA_m) .* RDimg;

photon_cts(photon_cts<0) = NaN;

% multiply with 
photon_noise_stdra = RDimg ./ sqrt(photon_cts);
photon_noise_stdif = photon_noise_stdra .* (pi ./ SFimg .* (d_au.^2));
photon_noise_mad_stdif = photon_noise_stdif .* norminv(0.75);

end

% photon_cts = TRRYRAdata.img * (A*sr*tsec) *6.55*10^(-3) ./ ((h_planck*c_light) ./ (WA.*10^(-9))) * eta;
% 
% std_rad = TRRYRAdata.img ./ sqrt(photon_cts);
% std_if = std_rad .* pi ./ SFdata.img .* (d_au.^2);
% 
% mad_photon_from_med = std_if .* norminv(0.75);
% 
% mad_photon_from_med = permute(mad_photon_from_med(:,:,bands),[3,1,2]);