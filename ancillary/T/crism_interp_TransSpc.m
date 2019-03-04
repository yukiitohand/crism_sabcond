function [T_interpd] = crism_interp_TransSpc(T,bands,RAdata)

[valid_samples] = crism_examine_valid_Columns(RAdata);

T_valid = T(bands,valid_samples);

T_interpd = nan(size(T));
for j=1:size(T_valid,1)
    %fprintf('j=%d\n',j);
    Tc = T_valid(j,:);
    [Tc_est] = naninterplocalLAD(Tc);
    T_interpd(bands(j),valid_samples) = Tc_est;
end