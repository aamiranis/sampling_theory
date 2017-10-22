%% load data

load('usps_all');
[k,l,m] = size(data); % k = 256, l = 1100, m = 10

rng('default');

for iter = 1:10
    
    fprintf('Generating data-set %d...\n', iter);
    
    % number of samples to consider from each class
    num_samples = 100;
    N = m*num_samples;

    % graph signal
    f = zeros(N,1);
    X = zeros(N,k);
    for i = 1:m
        f((i-1)*num_samples + 1 : (i-1)*num_samples + num_samples) = i;
        X((i-1)*num_samples + 1 : (i-1)*num_samples + num_samples, :) = data(:, randsample(1:l, num_samples), i)';
    end

    % membership functions
    mem_fn = false(N,m);
    for i = 1:m
        mem_fn(f==i,i) = true;
    end

    % compute pairwise distances
    distance = zeros(N);
    % Assign lower triangular part only in loop, saves time
    for i = 1:N
        for j = 1:i-1
            distance(i,j) = sqrt((X(i,:) - X(j,:))*(X(i,:) - X(j,:))');
        end
    end
    % Complete upper triangular part
    distance = distance + distance';

    % Number of nearest neighbors
    knn_param = 10;

    % Calculating distances of k-nearest neighbors
    knn_distance = zeros(N,1);
    for i = 1:N
        % sort all possible neighbors according to distance
        temp = sort(distance(i,:), 'ascend');
        % select k-th neighbor: knn_param+1, as the node itself is considered
        knn_distance(i) = temp(knn_param + 1);
    end

    % sparsification matrix
    nodes_to_retain = true(N);
    for i = 1:N
        nodes_to_retain(i, distance(i,:) > knn_distance(i) ) = false;
        nodes_to_retain(i,i) = false; % diagonal should be zero
    end

    nodes_to_retain( nodes_to_retain ~= nodes_to_retain' ) = true; 
    % do not symmetrize the knn matrix

    % computing sigma
    sigma = 1/3 * mean(knn_distance);

    % Creating adjacency matrix: only compute values for nodes to retain
    A = zeros(N);
    A(nodes_to_retain) = exp( -distance(nodes_to_retain).^2 / (2*sigma^2) );
%     A(nodes_to_retain) = 1;
    A = sparse(A);
    
    % find the largest SCC (because we need to define an irreducible random walk)
    [num_scc, scc_id] = graphconncomp(A);
    max_scc_sz = 0;
    max_scc_id = 0;
    for i = 1:num_scc
        scc_sz = sum(scc_id == i);
        if(max_scc_sz <= scc_sz)
            max_scc_sz = scc_sz;
            max_scc_id = i;
        end
    end
    scc = (scc_id == max_scc_id);
    A_scc = A(scc,scc);
    
    % save adjacency matrix
    save(['set', num2str(iter), '_dir.mat'],'A', 'scc', 'A_scc', 'mem_fn', 'X', 'sigma');

end
