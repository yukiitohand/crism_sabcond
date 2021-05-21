function [logmodel_edge_cor] = correct_edge_channels_bprep(wv,logmodel,logYif_isnan)


[B,L] = size(logmodel);
logmodel_edge_cor  = logmodel;
logYif_isgood = (logYif_isnan==0);
for l=1:L
%     b_lg = find(logYif_isgood(:,l),2,'first'); % left good
    b_rg = find(logYif_isgood(:,l),2,'last'); % right good
%     b_lb = 1:b_lg(1); % left bad (the most left good one is also considered)
    b_rb = b_rg(2):B; % right bad (the most right good oneis also considered)
%     b_lg4ex = b_lg(2):(b_lg(2)+1); % left good points for extrapolation.
    b_rg4ex = (b_rg(1)-1):b_rg(1); % right good points for extrapolation.
%     logmodel_edge_cor(b_lb,l) = interp1(wv(b_lg4ex),logmodel(b_lg4ex,l),wv(b_lb),'linear','extrap');
    logmodel_edge_cor(b_rb,l) = interp1(wv(b_rg4ex),logmodel(b_rg4ex,l),wv(b_rb),'linear','extrap');
end

end