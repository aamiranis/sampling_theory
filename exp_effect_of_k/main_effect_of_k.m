%% data and required toolboxes

addpath(genpath('../graphs'));
addpath(genpath('../sampling_methods'));
addpath(genpath('../reconstruction_methods'));

%% graph

fprintf('\ngraph = %d\n', graph)
fprintf('\nsignal = %d\n', signal)

rng('default');

switch graph
    case 1
        load erdos_renyi_graph_01 % contains A, Ln, deg, U, lambda
        conn_prob = 1;
        
    case 2
        load erdos_renyi_graph_05
        conn_prob = 5;
        
    case 3
        load erdos_renyi_graph_10
        conn_prob = 10;
end

N = length(A); % number of nodes

% higher power of Laplacian
k1 = 2;
Ln_k1 = Ln;
for i = 1:(k1-1)
    Ln_k1 = Ln_k1*Ln;
end
Ln_k1 = 0.5*(Ln_k1+Ln_k1.');

k2 = 8;
Ln_k2 = Ln;
for i = 1:(k2-1)
    Ln_k2 = Ln_k2*Ln;
end
Ln_k2 = 0.5*(Ln_k2+Ln_k2.');

k3 = 14;
Ln_k3 = Ln;
for i = 1:(k3-1)
    Ln_k3 = Ln_k3*Ln;
end
Ln_k3 = 0.5*(Ln_k3+Ln_k3.');


%% iterate over signals
num_sig = 50;

num_trials = 50; % average results for randomized methods

r = 50;

% number of nodes to be sampled
sample_size = 20:20:200;
num_points = length(sample_size);

error_set_L_k1 = zeros(num_points, num_sig);
error_set_L_k2 = zeros(num_points, num_sig);
error_set_L_k3 = zeros(num_points, num_sig);

samples_L_k1 = cell(num_points,1); 
samples_L_k2 = cell(num_points,1); 
samples_L_k3 = cell(num_points,1); 

prev_sample_size = 0;

prev_samples_L_k1 = false(N,1);
prev_samples_L_k2 = false(N,1);
prev_samples_L_k3 = false(N,1);

cond_L_k1 = zeros(num_points, 1);
cond_L_k2 = zeros(num_points, 1);
cond_L_k3 = zeros(num_points, 1);

UR = U(:,1:r);

for i = 1:num_points
    
    fprintf('******************\n\ni = %d\n\n', i);
    
    % number nodes to be sampled in the current iteration
    sample_size_inc = sample_size(i) - prev_sample_size;
    
    % sample and reconstruct
    
    %%%% max lamba_min[(L^k)_S^c] with k = k1
    samples_L_k1{i} = compute_S_L_k(Ln_k1, k1, sample_size_inc, prev_samples_L_k1);
    cond_L_k1(i) = cond(UR(samples_L_k1{i},:));
    fprintf('number of nodes sampled (with k = %d) = %d \n', k1, sum(samples_L_k1{i}));
    fprintf('condition number (with k = %d) = %f \n\n', k1, cond_L_k1(i));
    
    %%%% max lamba_min[(L^k)_S^c] with k = k2
    samples_L_k2{i} = compute_S_L_k(Ln_k2, k2, sample_size_inc, prev_samples_L_k2);
    cond_L_k2(i) = cond(UR(samples_L_k2{i},:));
    fprintf('number of nodes sampled (with k = %d) = %d \n', k2, sum(samples_L_k2{i}));
    fprintf('condition number (with k = %d) = %f \n\n', k2, cond_L_k2(i));
    
    %%%% max lamba_min[(L^k)_S^c] with k = k3
    samples_L_k3{i} = compute_S_L_k(Ln_k3, k3, sample_size_inc, prev_samples_L_k3);
    cond_L_k3(i) = cond(UR(samples_L_k3{i},:));
    fprintf('number of nodes sampled (with k = %d) = %d \n', k3, sum(samples_L_k3{i}));
    fprintf('condition number (with k = %d) = %f \n\n', k3, cond_L_k3(i));
    
    
    % update variable for next iteration
    prev_sample_size = sample_size(i);
    
    prev_samples_L_k1 = samples_L_k1{i};
    prev_samples_L_k2 = samples_L_k2{i};
    prev_samples_L_k3 = samples_L_k3{i};
    
end


for sig = 1:num_sig
%% generate signal

    switch signal
        case 1 % low noise signal
            % exact bandlimited
%             r = 60; % floor(0.1*N); % cutoff frequency index
            omega = lambda(r); % bandwidth
            c = 1+0.5*randn([r,1]); % low pass coefficients
            d = zeros(N-r,1); % exact bandlimited + noise
            f = U*[c; d];
            f = f/norm(f);
            fn = f;
            % add noise (error computation with noiseless signal!)
%             snr = 40;
%             sigma = norm(f)/(10^(snr/20));
%             fn = f + sigma*randn(N,1);
        case 2 % higher noise signal
            % exact bandlimited
%             r = 60; % floor(0.1*N); % cutoff frequency index
            omega = lambda(r); % bandwidth
            c = 1+0.5*randn([r,1]); % low pass coefficients
            d = zeros(N-r,1); % exact bandlimited + noise
            f = U*[c; d];
            f = f/norm(f);
            % add noise (error computation with noiseless signal!)
            snr = 20;
            sigma = norm(f)/(10^(snr/20));
            fn = f + sigma*randn(N,1);
        case 3
            % random signal with exponentially decaying gft coeffs after lambda_r;
%             r = 60;
            c = 1 + 0.5*randn(N,1);
            c(r:end) = exp(-4*(lambda(r:end) - lambda(r))).*c(r:end);
            f = U*c;
            f = f/norm(f);
            fn = f;
    end
    
    %% variables for looping and storing results

    error_L_k1 = zeros(num_points, 1);
    error_L_k2 = zeros(num_points, 1);
    error_L_k3 = zeros(num_points, 1);

    %% main loop: different reconstruction methods for different set selections
    
    for i = 1:num_points

        fprintf('******************\n\ni = %d\n\n', i);

        % sample and reconstruct

        %%%% k = k1 + exact BL reconstruction %%%%
        error_L_k1(i) = reconstruction_exact_BL(f, fn, samples_L_k1{i}, U(:,1:r));
        fprintf('reconstruction MSE (with k = %d) = %f \n\n', k1, error_L_k1(i));
        
        %%%% k = k2 + exact BL reconstruction %%%%
        error_L_k2(i) = reconstruction_exact_BL(f, fn, samples_L_k2{i}, U(:,1:r));
        fprintf('reconstruction MSE (with k = %d) = %f \n\n', k2, error_L_k2(i));

        %%%% k = k3 + exact BL reconstruction %%%%
        error_L_k3(i) = reconstruction_exact_BL(f, fn, samples_L_k3{i}, U(:,1:r));
        fprintf('reconstruction MSE (with k = %d) = %f \n\n', k3, error_L_k3(i));

    end

    error_set_L_k1(:,sig) = error_L_k1;
    error_set_L_k2(:,sig) = error_L_k2;
    error_set_L_k3(:,sig) = error_L_k3;
    
%     save(['results/test_result_G' num2str(graph) '_P' num2str(conn_prob) '_F' num2str(signal) '.mat'],'sample_size',...
%     'samples_L_k1', 'samples_L_k2', 'samples_L_k3', 'error_set_L_k1','error_set_L_k2','error_set_L_k3');

end
    
save(['results/result1_ER' '_P' num2str(conn_prob) '_F' num2str(signal) '.mat'],'sample_size',...
    'samples_L_k1', 'samples_L_k2', 'samples_L_k3', 'error_set_L_k1','error_set_L_k2','error_set_L_k3');

%% plots

% figure, shadedErrorBar(sample_size, mean(error_set_L_k,2), std(error_set_L_k,0,2), '-b');
% hold on 
% shadedErrorBar(sample_size, mean(error_set_U_SR,2), std(error_set_U_SR,0,2), '-r');
% shadedErrorBar(sample_size, mean(error_set_S_U,2), std(error_set_S_U,0,2), '-c');
% shadedErrorBar(sample_size, mean(error_set_rand_U_VR,2), std(error_set_rand_U_VR,0,2), '-k');
% shadedErrorbar(sample_size, mean(error_set_rand_uni,2), std(error_set_rand_uni,0,2),'-g');
% xlabel('Sample size')
% ylabel('Reconstruction MSE')
% legend('Proposed method', 'M1', 'M2', 'M3', 'random sampling')

% figure, plot(sample_size, mean(error_set_L_k,2), '-b');
% hold on 
% plot(sample_size, mean(error_set_U_SR,2), '-r');
% plot(sample_size, mean(error_set_S_U,2), '-c');
% plot(sample_size, mean(error_set_rand_U_VR,2), '-k');
% % plot(sample_size, mean(error_set_rand_uni,2), '-g');
% xlabel('Sample size')
% ylabel('Reconstruction MSE')
% legend('Proposed method', 'M1', 'M2', 'M3', 'random sampling')
% title('G1 with F2')

