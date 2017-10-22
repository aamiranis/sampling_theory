%% data and required toolboxes

addpath(genpath('../graphs'));
addpath(genpath('../sampling_methods'));
addpath(genpath('../reconstruction_methods'));

%% graph

% number of nodes to be sampled
sample_size = 60:10:150;
num_points = length(sample_size);

error_L_k = zeros(num_points, 10);
error_U_SR = zeros(num_points, 10);
error_S_U = zeros(num_points, 10);
error_rand_uni = zeros(num_points, 10);

for set = 1:10

    load(['set' num2str(set) '_dir.mat']);

    A = A(scc,scc);
    X = X(scc,:);
    mem_fn = mem_fn(scc,:);
    f = mem_fn;
    fn = f;

    % try normalizing the adjacency
    d = sum(A,2);
    A = diag(1./d)*A;

    num_classes = size(mem_fn,2);
    N = length(A); % number of nodes

    % EVD of A
    [V,mu] = eig(full(A));
    mu = diag(mu);
    mu_max = max(abs(mu));

    % variation operator
    Ln = speye(N) - A/mu_max;
    % EVD of Ln
    lambda = 1 - mu/mu_max;
    % [V,lambda] = eig(full(Ln));
    % lambda = diag(lambda);
    [~,id] = sort(abs(lambda),'ascend');
    lambda = lambda(id);

    V = V(:,id); % inverse GFT mtx
    U = inv(V); % GFT mtx

    % higher power of Laplacian
    k = 4;
    Ln_k = Ln;
    for i = 1:(k-1)
        Ln_k = Ln_k*Ln;
    end
    % Ln_k_sym = 0.5*(Ln_k+Ln_k.');
    Ln_k_sym = Ln_k'*Ln_k;

    %% iterate over signals

    num_trials = 100;
    
    samples_L_k = cell(num_points,1); 
    samples_U_SR = cell(num_points,1); 
    samples_S_U = cell(num_points,1);

    prev_sample_size = 0;

    prev_samples_L_k = false(N,1);
    prev_samples_S_U = false(N,1);
    prev_samples_U_SR = false(N,1);

    cond_L_k = zeros(num_points, 1);
    cond_U_SR = zeros(num_points, 1);
    cond_S_U = zeros(num_points, 1);

    gft = U*mem_fn;
    % energy_fraction = sqrt(sum(cumsum(abs(gft).^2),2))/norm(gft,'fro');
    energy_fraction = sum(cumsum(abs(gft).^2),2) / (norm(gft,'fro')^2);
    r_test = find(energy_fraction>0.8, 1);
    % omega = lambda(r);

    r = 50;

    for i = 1:num_points

        fprintf('******************\n\ni = %d\n\n', i);

    %     r = sample_size(i);

        % number nodes to be sampled in the current iteration
        sample_size_inc = sample_size(i) - prev_sample_size;

        % sample and reconstruct

        %%%% max lamba_min[(L^k)_S^c]
        [samples_L_k{i}, omega, ~] = compute_S_L_k_dir_1(Ln_k_sym, 2*k, sample_size_inc, prev_samples_L_k);
        fprintf('number of nodes sampled (proposed) = %d \n\n', sum(samples_L_k{i}));
    %     error_L_k(i) = reconstruction_pocs_dir(f, fn, samples_L_k{i}, Ln, omega);
        error_L_k(i,set) = reconstruction_exact_BL_dir(f, fn, samples_L_k{i}, V(:,1:r));
        cond_L_k(i) = cond(V(samples_L_k{i},1:r));
        fprintf('condition number (proposed) = %f \n', cond_L_k(i));
        fprintf('reconstruction MSE (max lamba_min[(L^k)_S^c]) = %f \n\n', error_L_k(i,set));

        %%%% max sigma_min[V_SR]
        samples_U_SR{i} = compute_S_U_SR_dir(V(:,1:r), sample_size_inc, prev_samples_U_SR);
        fprintf('number of nodes sampled (M1) = %d \n\n', sum(samples_U_SR{i}));
        [error_U_SR(i,set),f_recon] = reconstruction_exact_BL_dir(f, fn, samples_U_SR{i}, V(:,1:r));
        cond_U_SR(i) = cond(V(samples_U_SR{i},1:r));
        fprintf('condition number (M1) = %f \n', cond_U_SR(i));
        fprintf('reconstruction MSE (max sigma_min[U_SR]) = %f \n\n', error_U_SR(i,set));

        %%%% max S^T V_VR (Shomorony)
        samples_S_U{i} = compute_S_S_U_dir(V, sample_size_inc, prev_samples_S_U);
        fprintf('number of nodes sampled (M2) = %d \n\n', sum(samples_S_U{i}));
        error_S_U(i,set) = reconstruction_exact_BL_dir(f, fn, samples_S_U{i}, V(:,1:r));
        cond_S_U(i) = cond(V(samples_S_U{i},1:r));
        fprintf('condition number (M2) = %f \n', cond_S_U(i));
        fprintf('reconstruction MSE (max S^T U_VR) = %f \n\n', error_S_U(i,set));

         %%%% randomized uniform sampling and exact BL reconstruction %%%%
        error = zeros(num_trials,1);
        for j = 1:num_trials
            error(j) = rand_sample_reconstruct_uni_dir(f, fn, V(:,1:r), sample_size(i));
        end
        error_rand_uni(i,set) = mean(error);
        fprintf('reconstruction MSE (randomized uniform) = %f \n\n', error_rand_uni(i,set));


        % update variable for next iteration
        prev_sample_size = sample_size(i);

        prev_samples_L_k = samples_L_k{i};
        prev_samples_S_U = samples_S_U{i};
        prev_samples_U_SR = samples_U_SR{i};
        
        save('results/usps_results.mat','sample_size','error_L_k', 'error_U_SR', 'error_S_U', 'error_rand_uni');
    end


    % save('results/usps_set_1_dir_results_norm.mat', 'error_L_k', 'error_U_SR', 'error_S_U', 'error_rand_uni',...
    %     'cond_L_k', 'cond_U_SR', 'cond_S_U',...
    %     'samples_L_k', 'samples_U_SR', 'samples_S_U');

    % save('results/usps_set_1_dir_results_est_cutoff.mat','error_L_k', 'samples_L_k');
    
    save('results/usps_results.mat','sample_size','error_L_k', 'error_U_SR', 'error_S_U', 'error_rand_uni');

end

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

