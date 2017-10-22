function [ S_opt, count ] = compute_S_L_k_lobpcg( L, prec_fun, k, num_nodes_to_add, current_S_opt )
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

omega = 0;

function result = operatorA(x)
    M = size(x,2);
    result = zeros(N,M);
    
    if (M == 1)
        result(~S_opt) = x;
    else
        result(~S_opt,:) = x;
    end
    
    for i = 1:k
%         result = (L * result + result) / 3;
        result = L * result;
    end
    
    if (M == 1)
        result = result(~S_opt);
    else
        result = result(~S_opt,:);
    end
    
%     result = result - omega * x;
    count = count + k;
end

function x = operatorT(x)
    precon_tol = 1e-4;
    max_iter_precon = N/10;
    
    k_cg = 0;

    r_cg = x - operatorA(x);
    p_cg = r_cg;
    rr_cg_old = r_cg'*r_cg;

    while ( (sqrt(rr_cg_old) > precon_tol) && (k_cg <= max_iter_precon) )
        k_cg = k_cg + 1;

        A_times_p_int = operatorA(p_cg);

        alpha_cg = rr_cg_old / (p_cg' * A_times_p_int);
        x = x + alpha_cg * p_cg;
        r_cg = r_cg - alpha_cg * A_times_p_int;

        rr_cg_new = r_cg'*r_cg;
        p_cg = r_cg + (rr_cg_new/rr_cg_old) * p_cg;

        rr_cg_old = rr_cg_new;
    end
end

% function result = operatorT(x)
%     M = size(x,2);
%     result = zeros(N,M);
%     
%     if (M == 1)
%         result(~S_opt) = x;
%     else
%         result(~S_opt,:) = x;
%     end
%     
%     for i = 1:k
%         result = prec_fun(result);
%     end
%     
%     if (M == 1)
%         result = result(~S_opt);
%     else
%         result = result(~S_opt,:);
%     end
% end


% iterations_for_convergence = zeros(num_nodes_to_add,1);

for iter = 1:num_nodes_to_add

%     tic

%   fprintf('Iteration %d of %d...\n', iter, num_nodes_to_add);

    % create index vector for Sc from indicator functions
    q = find(~S_opt);

    % compute minimum eigen-pair: efficient way
    initial_x = ones(length(q),1);
    initial_x = initial_x / norm(initial_x);
    
    [y, ~, failure_flag] = lobpcg(initial_x, @(x)operatorA(x), 1e-3, N);
%     [y, omega, failure_flag, lambda_history, residual_norm_history] = lobpcg(initial_x, @(x)operatorA(x), [], @(x)operatorT(x), 1e-4, N);
    
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

