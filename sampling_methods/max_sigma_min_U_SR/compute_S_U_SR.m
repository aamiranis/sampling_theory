function S = compute_S_U_SR(UR, sample_size_inc, prev_S)


S = prev_S;

for i = 1:sample_size_inc
    
    unsampled_nodes = find(~prev_S);
    max_sigma_min = 0;
    
    for j = 1:length(unsampled_nodes)
        
        temp_S = prev_S;
        temp_S(unsampled_nodes(j)) = 1;

        new_sigma_min = min(svd(UR(temp_S,:)));

        if(new_sigma_min >= max_sigma_min)
            max_sigma_min = new_sigma_min;
            S = temp_S;
        end
        
    end
    
    % update vars before choosing the next node
    prev_S = S;
   
end

