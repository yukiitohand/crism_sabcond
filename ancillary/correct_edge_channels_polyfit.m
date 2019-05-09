function [logmodel_edge_cor] = correct_edge_channels_polyfit(wv,logmodel)

nl=2;
% nr=3;

[B,L] = size(logmodel);
logmodel_edge_cor  = logmodel;

b_lb = 1:7; % left bad (the most left good one is also considered)
% b_rb = (B-4):B; % right bad (the most right good oneis also considered)
b_lg4ex = 8:20; % left good points for extrapolation.
% b_rg4ex = (B-19):(B-5); % right good points for extrapolation.

wv_micron = wv/1000;
wv_mat_l = (wv_micron).^(nl:-1:0);
% wv_mat_r = (wv_micron).^(nr:-1:0);

for l=1:L
    logmodel_edge_param_l = polyfit(wv_micron(b_lg4ex),logmodel(b_lg4ex,l),nl);
%     logmodel_edge_param_r = polyfit(wv_micron(b_rg4ex),logmodel(b_rg4ex,l),nr);
    logmodel_edge_cor(b_lb,l) = wv_mat_l(b_lb,:) * logmodel_edge_param_l';
%     logmodel_edge_cor(b_rb,l) = wv_mat_r(b_rb,:) * logmodel_edge_param_r';
end

end