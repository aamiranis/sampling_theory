%% data and required toolboxes

addpath(genpath('../graphs'));
addpath(genpath('../sampling_methods'));
addpath(genpath('../lobpcg'));

%% graph

rng('default');

switch graph
    case 1
        load erdos_renyi_graph_1k % contains A, Ln, deg, U, lambda
    case 2
        load erdos_renyi_graph_5k
    case 3
        load erdos_renyi_graph_10k
    case 4
        load erdos_renyi_graph_20k
    case 5
        load erdos_renyi_graph_50k
    case 6
%         load erdos_renyi_graph_100k
        erdos_renyi_diff_sizes
end

N = length(A); % number of nodes

%% iterate over signals

num_reps = 1;

% number of nodes to be sampled
% sample_size_list = [10 20 30 40 50];
% sample_size_list = 50;
% sample_size_list = floor([0.2 0.4 0.6 0.8 1.0]/100 * N);
% sample_size_list = floor([0.5 1 1.5 2 2.5]/100 * N);
sample_size_list = floor([1 2 3 4 5]/100 * N);

samples_S_U = false(N,1);
samples_U_SR = false(N,1);
samples_L_4 = false(N,1);
samples_L_6 = false(N,1);
samples_L_8 = false(N,1);

time_U_SR = zeros(length(sample_size_list),1);
time_S_U = zeros(length(sample_size_list),1);
time_L_4 = zeros(length(sample_size_list),1);
time_L_6 = zeros(length(sample_size_list),1);
time_L_8 = zeros(length(sample_size_list),1);

tic  
L_prec = ichol(Ln); prec_fun = @(x)L_prec'\(L_prec\x);
initial_x = randn(N, sample_size_list(end));
for j = 1:size(initial_x,2)
    initial_x(:,j) = initial_x(:,j)/norm(initial_x(:,j));
end
[UR, ~, failure_flag, ~, residual_norm_history] = lobpcg(initial_x, Ln, 1e-3, N);
% [UR, ~, failure_flag, ~, residual_norm_history] = lobpcg(initial_x, Ln, [], @(x)prec_fun(x), 1e-4, N);
if (failure_flag == 1)
    fprintf('Did not converge!...\n');
end
t0 = toc;

prev_sample_size = 0;

for i = 1:length(sample_size_list)
    
    sample_size_inc = sample_size_list(i) - prev_sample_size;
    prev_sample_size = sample_size_list(i);
    
    tic
    samples_S_U = compute_S_S_U(UR, sample_size_inc, samples_S_U);
    if (i == 1)
        time_S_U(i) = t0 + toc;
    else
        time_S_U(i) = time_S_U(i-1) + toc;
    end
    
%     tic
%     samples_U_SR = compute_S_U_SR(UR, sample_size_inc, samples_U_SR);
%     if (i == 1)
%         time_U_SR(i) = t0 + toc;
%     else
%         time_U_SR(i) = time_U_SR(i-1) + toc;
%     end
    
    %%%% max lamba_min[(L^k)_S^c]
    k = 4;
    tic
    [samples_L_4, count_4] = compute_S_L_k_lobpcg(Ln, prec_fun, k, sample_size_inc, samples_L_4);
%     [samples_L_4, count_4] = compute_S_L_k_lobpcg_proj(Ln, prec_fun, k, sample_size_inc, samples_L_4);
    if (i == 1)
        time_L_4(i) = toc;
    else
        time_L_4(i) = time_L_4(i-1) + toc;
    end
    
    k = 6;
    tic
    [samples_L_6, count_6] = compute_S_L_k_lobpcg(Ln, prec_fun, k, sample_size_inc, samples_L_6);
%     [samples_L_6, count_6] = compute_S_L_k_lobpcg_proj(Ln, prec_fun, k, sample_size_inc, samples_L_6);
    if (i == 1)
        time_L_6(i) = toc;
    else
        time_L_6(i) = time_L_6(i-1) + toc;
    end
    

    k = 8;
    tic
    [samples_L_8, count_8] = compute_S_L_k_lobpcg(Ln, prec_fun, k, sample_size_inc, samples_L_8);
%     [samples_L_8, count_8] = compute_S_L_k_lobpcg_proj(Ln, prec_fun, k, sample_size_inc, samples_L_8);
    if (i == 1)
        time_L_8(i) = toc;
    else
        time_L_8(i) = time_L_8(i-1) + toc;
    end
    
    fprintf('Samples added = %d\n', sample_size_list(i));
end