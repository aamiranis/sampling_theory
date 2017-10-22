function S = compute_S_S_U(U, sample_size_inc, prev_S)

m = sum(prev_S); % size of previous sampling set.
I = speye(length(prev_S));

S = prev_S;

% unsampled_nodes = find(~prev_S);

for i = 1:sample_size_inc
    % represent U(:,m+1) as a combination of U(:,1:m) and eye(:,~prev_S)
    B = sparse([U(:, 1:m), I(:,~prev_S)]);
    c = B \ U(:,m+1);
    [~, id] = max(abs(c(m+1:end)));

    unsampled_nodes = find(~prev_S);
    v = unsampled_nodes(id); % node to sample
    S(v) = 1;
    
    m = m+1;
    prev_S = S;
end


