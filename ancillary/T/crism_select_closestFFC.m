function [ dirname_ffc,dirname_ffc_mcb,dirname_ffc_mca ] = crism_select_closestFFC(sclk)
% [ dirname_ffc,dirname_ffc_mcb,dirname_ffc_mca ] = crism_select_closestFFC(sclk)
%   obtain the dirname of the closest FFC measurements, accompanied with
%   the most recent before and most closest post FFC.
%  INPUT
%    sclk: double, sclk time
%  OUTPUTS
%    dirname_ffc: string, dirname of the closest FFC
%    dirname_ffc_mcb: string, dirname of the most closest prior FFC
%    dirname_ffc_mca: string, dirname of the most closest post FFC
%
%    dirnames are in the form of
%       cccnnnnnn
%        ccc = class type of the obervation
%        nnnnnnn = observation id 
%  

load('ffc_info.mat','ffc_info');

sclkFFC01_starts = [ffc_info.sclk_start_01];
sclkFFC01_stops  = [ffc_info.sclk_stop_01];

sclkFFC01 = (sclkFFC01_starts+sclkFFC01_stops)/2;

t_diff = sclk-sclkFFC01;
t_diff_abs = abs(t_diff);

[~,idx_closest] = nanmin(t_diff_abs);


if t_diff(idx_closest) < 0
    idx_mcb = idx_closest-1;
    idx_mca = idx_closest;
else
    idx_mcb = idx_closest;
    idx_mca = idx_closest+1;
end

if idx_mcb==0
    idx_mcb = [];
end
if idx_mca > length(sclkFFC01)
    idx_mca = [];
end

dirname_ffc = ffc_info(idx_closest).dirname;

if isempty(idx_mcb)
    dirname_ffc_mcb = '';
else
    dirname_ffc_mcb = ffc_info(idx_mcb).dirname;
end

if isempty(idx_mca)
    dirname_ffc_mca = '';
else
    dirname_ffc_mca = ffc_info(idx_mca).dirname;
end



end