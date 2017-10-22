function [ S_opt, count ] = compute_S_L_k_lobpcg_proj( L, prec_fun, k, num_nodes_to_add, current_S_opt )
%   AUTHOR: Aamir Anis, USC
%   This function computes the optimal sampling set of a given size 
%   "S_opt_size" that maximizes the cutoff frequency.

% % %
% PARAMETER DESCRIPTION
% 
% INPUT
% L_k: kth power of Laplacian
% S_opt_size:  Desired size of the optimal set
% k: Power of Laplacian while computing cutoff, higher the order,
% greater the accuracy, but the complexity is also higher.
% 
% OUTPUT
% S_opt: Optimal set as a logical vector
% omega: cutoff of the optimal set
% omega_list: List of computed cutoffs
% 
% % %

% fprintf('Starting optimal set search...\n');
N = length(L);

count = 0;

% Initialization : If previous state available, initialize to that
if (exist('current_S_opt','var'))
    S_opt = current_S_opt;
else
    S_opt = false(N,1);
end

function x = operatorA(x)
    for i = 1:k
        x = L * x;
    end
    count = count + k;
end

function x = operatorT(x)
    for i = 1:k
        x = prec_fun(x);
    end
end

% iterations_for_convergence = zeros(num_nodes_to_add,1);

for iter = 1:num_nodes_to_add

%     tic

%   fprintf('Iteration %d of %d...\n', iter, num_nodes_to_add);

    % create index vector for Sc from indicator functions
    q = find(~S_opt);

    % compute minimum eigen-pair: efficient way
    initial_x = ones(length(q),1);
    initial_x = initial_x / norm(initial_x);
    [ y, omega, failure_flag, log ] = eig_lopcg_proj( initial_x, ~S_opt, @(x)operatorA(x), @(x)operatorT(x), 1e-4, 500);
%     log(:,3)
    failure_flag = abs(1 - failure_flag); 
    
    if (failure_flag == 1)
        fprintf('k = %d did not converge while adding node %d.\n', k, iter);
    end
        
    % find direction of maximum increment in reduced (|Sc|) dimensions
    [~,max_index] = max(abs(y));

    % Find corresponding node in N dimensions
    node_to_add = q(max_index);

    % Update indicator function
    S_opt(node_to_add) = 1;

%     fprintf('Nodes added = %d...\n', sum(S_opt));
    
%     toc

end

% fprintf('Finished.\n');



end

