%% data and required toolboxes

addpath(genpath('../graphs'));
addpath(genpath('../sampling_methods'));
addpath(genpath('../reconstruction_methods'));

%% graph

rng('default');

switch graph
    case 1
        load erdos_renyi_graph_1 % contains A, Ln, deg, U, lambda
    case 2
        load small_world_graph
    case 3
        load barbasi_albert_graph
end

N = length(A); % number of nodes

% higher powers of Laplacian
k = 2;
Ln_2 = Ln;
for i = 1:(k-1)
    Ln_2 = Ln_2*Ln;
end
Ln_2 = 0.5*(Ln_2+Ln_2.');

k = 8;
Ln_8 = Ln;
for i = 1:(k-1)
    Ln_8 = Ln_8*Ln;
end
Ln_8 = 0.5*(Ln_8+Ln_8.');

k = 14;
Ln_14 = Ln;
for i = 1:(k-1)
    Ln_14 = Ln_14*Ln;
end
Ln_14 = 0.5*(Ln_14+Ln_14.');

%% iterate over signals
num_sig = 50;

num_trials = 50; % average results for randomized methods

r = 50;

% number of nodes to be sampled
sample_size = 20:20:200;
num_points = length(sample_size);

error_set_L_2 = zeros(num_points, num_sig);
error_set_L_8 = zeros(num_points, num_sig);
error_set_L_14 = zeros(num_points, num_sig);
error_set_U_SR = zeros(num_points, num_sig);
error_set_S_U = zeros(num_points, num_sig);

error_set_rand_U_VR = zeros(num_points, num_sig); 
error_set_rand_uni = zeros(num_points, num_sig);

samples_L_2 = cell(num_points,1);
samples_L_8 = cell(num_points,1);
samples_L_14 = cell(num_points,1);
samples_U_SR = cell(num_points,1); 
samples_S_U = cell(num_points,1);

prev_sample_size = 0;

prev_samples_L_2 = false(N,1);
prev_samples_L_8 = false(N,1);
prev_samples_L_14 = false(N,1);
prev_samples_U_SR = false(N,1);
prev_samples_S_U = false(N,1);

cond_L_2 = zeros(num_points, 1);
cond_L_8 = zeros(num_points, 1);
cond_L_14 = zeros(num_points, 1);
cond_U_SR = zeros(num_points, 1);
cond_S_U = zeros(num_points, 1);

UR = U(:,1:r);

for i = 1:num_points
    
    fprintf('******************\n\ni = %d\n\n', i);
    
    % number nodes to be sampled in the current iteration
    sample_size_inc = sample_size(i) - prev_sample_size;
    
    % sample and reconstruct
    
    %%%% max lamba_min[(L^k)_S^c]
    samples_L_2{i} = compute_S_L_k(Ln_2, 2, sample_size_inc, prev_samples_L_2);
    cond_L_2(i) = cond(UR(samples_L_2{i},:));
    fprintf('number of nodes sampled (proposed) = %d \n', sum(samples_L_2{i}));
    fprintf('condition number (proposed) = %f \n\n', cond_L_2(i));
    
    %%%% max lamba_min[(L^k)_S^c]
    samples_L_8{i} = compute_S_L_k(Ln_8, 8, sample_size_inc, prev_samples_L_8);
    cond_L_8(i) = cond(UR(samples_L_8{i},:));
    fprintf('number of nodes sampled (proposed) = %d \n', sum(samples_L_8{i}));
    fprintf('condition number (proposed) = %f \n\n', cond_L_8(i));
    
    %%%% max lamba_min[(L^k)_S^c]
    samples_L_14{i} = compute_S_L_k(Ln_14, 14, sample_size_inc, prev_samples_L_14);
    cond_L_14(i) = cond(UR(samples_L_14{i},:));
    fprintf('number of nodes sampled (proposed) = %d \n', sum(samples_L_14{i}));
    fprintf('condition number (proposed) = %f \n\n', cond_L_14(i));
    
    %%%% max sigma_min[U_SR]
    samples_U_SR{i} = compute_S_U_SR(U(:,1:r), sample_size_inc, prev_samples_U_SR);
    cond_U_SR(i) = cond(UR(samples_U_SR{i},:));
    fprintf('number of nodes sampled (M1) = %d \n', sum(samples_U_SR{i}));
    fprintf('condition number (M1) = %f \n\n', cond_U_SR(i));
    
    %%%% max S^T U_VR
    samples_S_U{i} = compute_S_S_U(U, sample_size_inc, prev_samples_S_U);
    cond_S_U(i) = cond(UR(samples_S_U{i},:));
    fprintf('number of nodes sampled (M2) = %d \n', sum(samples_S_U{i}));
    fprintf('condition number (M2) = %f \n\n', cond_S_U(i));
    
    % update variable for next iteration
    prev_sample_size = sample_size(i);
    
    prev_samples_L_2 = samples_L_2{i};
    prev_samples_L_8 = samples_L_8{i};
    prev_samples_L_14 = samples_L_14{i};
    prev_samples_U_SR = samples_U_SR{i};
    prev_samples_S_U = samples_S_U{i};
    
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

    error_L_2 = zeros(num_points, 1);
    error_L_8 = zeros(num_points, 1);
    error_L_14 = zeros(num_points, 1);
    error_U_SR = zeros(num_points, 1);
    error_S_U = zeros(num_points, 1);

    error_rand_U_VR = zeros(num_points, 1); 
    error_rand_uni = zeros(num_points,1);

    %% main loop: different reconstruction methods for different set selections
    
    for i = 1:num_points

        fprintf('******************\n\ni = %d\n\n', i);

        % sample and reconstruct

        %%%% max lamba_min[(L^k)_S^c] + pocs reconstruction %%%%
%         error_L_k(i) = reconstruction_pocs(f, fn, samples_L_k{i}, Ln, omega);
        error_L_2(i) = reconstruction_exact_BL(f, fn, samples_L_2{i}, U(:,1:r));
        fprintf('reconstruction MSE (max lamba_min[(L^k)_S^c]) = %f \n\n', error_L_2(i));
        
        %%%% max lamba_min[(L^k)_S^c] + pocs reconstruction %%%%
%         error_L_k(i) = reconstruction_pocs(f, fn, samples_L_k{i}, Ln, omega);
        error_L_8(i) = reconstruction_exact_BL(f, fn, samples_L_8{i}, U(:,1:r));
        fprintf('reconstruction MSE (max lamba_min[(L^k)_S^c]) = %f \n\n', error_L_8(i));

        %%%% max lamba_min[(L^k)_S^c] + pocs reconstruction %%%%
%         error_L_k(i) = reconstruction_pocs(f, fn, samples_L_k{i}, Ln, omega);
        error_L_14(i) = reconstruction_exact_BL(f, fn, samples_L_14{i}, U(:,1:r));
        fprintf('reconstruction MSE (max lamba_min[(L^k)_S^c]) = %f \n\n', error_L_14(i));
        
        %%%% max sigma_min[U_SR] + exact BL reconstruction (Chen et al.) %%%%
        [error_U_SR(i),f_recon] = reconstruction_exact_BL(f, fn, samples_U_SR{i}, U(:,1:r));
        fprintf('reconstruction MSE (max sigma_min[U_SR]) = %f \n\n', error_U_SR(i));
        
        %%%% max S^T U_VR + exact BL reconstruction (Shomorony et al.) %%%%
        error_S_U(i) = reconstruction_exact_BL(f, fn, samples_S_U{i}, U(:,1:r));
        fprintf('reconstruction MSE (max S^T U_VR) = %f \n\n', error_S_U(i));

        %%%% randomized U_VR based sampling and reconstruction (Chen et al.) %%%%
        error = zeros(num_trials,1);
        for j = 1:num_trials
%             error(j) = rand_sample_reconstruct_U_VR(f, fn, U(:,1:r), sample_size(i));
            error(j) = rand_sample_reconstruct_U_VR(f, fn, U(:,1:r), sample_size(i));
        end
        error_rand_U_VR(i) = mean(error);
        fprintf('reconstruction MSE (randomized U_VR based) = %f \n\n', error_rand_U_VR(i,1));

        %%%% randomized uniform sampling and exact BL reconstruction %%%%
        error = zeros(num_trials,1);
        for j = 1:num_trials
            error(j) = rand_sample_reconstruct_uni(f, fn, U(:,1:r), sample_size(i));
        end
        error_rand_uni(i) = mean(error);
        fprintf('reconstruction MSE (randomized uniform) = %f \n\n', error_rand_uni(i,1));

    end

    error_set_L_2(:,sig) = error_L_2;
    error_set_L_8(:,sig) = error_L_8;
    error_set_L_14(:,sig) = error_L_14;
    error_set_U_SR(:,sig) = error_U_SR;
    error_set_S_U(:,sig) = error_S_U;

    error_set_rand_U_VR(:,sig) = error_rand_U_VR; 
    error_set_rand_uni(:,sig) = error_rand_uni;


    save(['results/result_G' num2str(graph) '_F' num2str(signal) '.mat'],'sample_size','error_set_L_2',...
        'error_set_L_8','error_set_L_14','error_set_U_SR','error_set_S_U','error_set_rand_U_VR','error_set_rand_uni');

end
    
save(['results/result_G' num2str(graph) '_F' num2str(signal) '.mat'],'sample_size','error_set_L_2',...
    'error_set_L_8','error_set_L_14','error_set_U_SR','error_set_S_U','error_set_rand_U_VR','error_set_rand_uni');
