function [x_width] = wavelength_sample_get_width(x)
% [x_width] = wavelength_sample_get_width(x)
%  compute interval for x
x = x(:);
x_extend = zeros([length(x)+2,1]);
x_extend(2:end-1) = x;
x_extend(1) = 2*x(1)-x(2);
x_extend(end)=2*x(end)-x(end-1);
x_between = (x_extend(2:end) + x_extend(1:end-1))/2;
x_width = x_between(:,:,2:end) - x_between(:,:,1:end-1);

