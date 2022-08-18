function [x_est] = naninterplocalLAD(x)

isnan_x = isnan(x);
bad_idx = find(isnan_x);

x_est = x;

for bi=bad_idx
    %if bi==1
    %    fprintf('bi=%d\n',bi);
    %end
    hw = 4;
    is_suf_neighbor = false;
    while ~is_suf_neighbor
        left_idx = max(1,bi-hw);
        right_idx = min(bi+hw,length(x));
        x_local_idx = left_idx:right_idx;
        x_local = x(x_local_idx);
    
        x_local_good_bool = ~isnan(x_local);
        if sum(x_local_good_bool)<2
            hw = hw+1;
        else
            is_suf_neighbor=true;
        end
    end  
        
    x_local_good_idx = x_local_idx(x_local_good_bool);
    x_local_good = x_local(x_local_good_bool);
    
    A = [x_local_good_idx(:) ones(length(x_local_good_idx),1)];
    
    coef = lad_gadmm_a_v2(A,x_local_good(:));
    
    x_est(bi) = coef(1)*bi+coef(2);
end