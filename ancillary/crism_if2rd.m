function [RD] = crism_if2rd(IoF,SF,lbl)
% [RD] = crism_if2rd(IoF,SFimg,r)
%   convert I/F to radiance
%  Input Parameters
%   IoF: I/F 
%   SFimg: CDR SF
%   lbl : CDR LBL struct having "SOLAR_DISTANCE" field
%   Note: IoF and SFimg need to have sizes compatible with implicit expansion.
%  Output Parameters
%   RD : Radiance cube [L x S x B]

d_km = lbl.SOLAR_DISTANCE.value;
[ d_au ] = km2au( d_km );
RD = IoF .* ( SF / (pi*(d_au)^2) );

end