function [X] = interp_nan_column(X,Xisnan,wv)
L = size(X,1);
rt_nan = sum(Xisnan,1)/L;
for l=1:size(X,2)
    if any(Xisnan(:,l)) && rt_nan(l)<0.8
        good_idx = ~Xisnan(:,l);
        X(Xisnan(:,l),l) = interp1(wv(good_idx),X(good_idx,l),...
            wv(Xisnan(:,l)),'linear','extrap');
    end
end