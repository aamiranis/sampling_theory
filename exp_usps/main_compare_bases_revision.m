clear
%% data and required toolboxes

addpath(genpath('../graphs'));
addpath(genpath('../sampling_methods'));
addpath(genpath('../reconstruction_methods'));

%% graph

% number of nodes to be sampled
sample_size = 60:10:150;
num_points = length(sample_size);

error_adj = zeros(num_points, 1);
error_ha = zeros(num_points, 1);
error_rw = zeros(num_points, 1);

for set = 1:10

    load(['set' num2str(set) '_dir.mat']);

    A = A(scc,scc);
    X = X(scc,:);
    mem_fn = mem_fn(scc,:);
    f = mem_fn;
    fn = f;

    % try normalizing the adjacency
    d = sum(A,2);
    A_norm = diag(1./d)*A;

    num_classes = size(mem_fn,2);
    N = length(A); % number of nodes

    % EVD of A
    [V,mu] = eig(full(A_norm));
    mu = diag(mu);
    mu_max = max(abs(mu));

    %% adjancency based variation operator

    % adj based variation operator
    Ln_adj = speye(N) - A_norm/mu_max; % A/mu_max;

    % EVD of Ln
    lambda = 1 - mu/mu_max;
    [~,id] = sort(abs(lambda),'ascend');
    lambda_adj = lambda(id);
    V_adj = V(:,id); % inverse GFT mtx
    U_adj = inv(V_adj); % GFT mtx

    % higher power of Ln_adj
    k = 4;
    Ln_adj_k = Ln_adj;
    for i = 1:(k-1)
        Ln_adj_k = Ln_adj_k*Ln_adj;
    end
    Ln_adj_k_sym = Ln_adj_k'*Ln_adj_k;

    %% hub-authority based variation operator

    d_q = sum(A,2); % out deg
    d_q_inv = zeros(N,1);
    d_q_inv(d_q ~= 0) = 1./sqrt(d_q);
    D_q_inv = diag(d_q_inv);

    d_p = sum(A,1);
    d_p_inv = zeros(N,1);
    d_p_inv(d_p ~= 0) = 1./sqrt(d_p);
    D_p_inv = diag(d_p_inv);

    T = D_q_inv * A * D_p_inv;
    Ln_a = eye(N) - T'*T;
    Ln_h = eye(N) - T*T';
    Ln_ha = 0.5*(Ln_a + Ln_h);
    Ln_ha = 0.5*(Ln_ha + Ln_ha'); % ensuring symmetry

    [V_ha,lambda_ha] = eig(full(Ln_ha));
    [lambda_ha,id] = sort(diag(lambda_ha),'ascend');
    V_ha = V_ha(:,id);
    U_ha = V_ha';

    % higher power of Ln_ha
    k = 8;
    Ln_ha_k = Ln_ha;
    for i = 1:(k-1)
        Ln_ha_k = Ln_ha_k*Ln_ha;
    end
    Ln_ha_k = 0.5*(Ln_ha_k+Ln_ha_k.');

    %% random walk based variation operator

    P = diag(1./d_q)*A; 
    [st_dist, eig_val] = eigs(P',1);

    D = diag(sqrt(st_dist));
    D_inv = zeros(N,1);
    D_inv(st_dist~=0) = 1./sqrt(st_dist); 
    D_inv = diag(D_inv);

    Wn = D*P*D_inv;
    Ln_rw = eye(N) - 0.5*(Wn + Wn');

    [V_rw,lambda_rw] = eig(full(Ln_rw));
    [lambda_rw,id] = sort(diag(lambda_rw),'ascend');
    V_rw = V_rw(:,id);
    U_rw = V_rw';

    % higher power of Ln_rw
    k = 8;
    Ln_rw_k = Ln_rw;
    for i = 1:(k-1)
        Ln_rw_k = Ln_rw_k*Ln_rw;
    end
    Ln_rw_k = 0.5*(Ln_rw_k+Ln_rw_k.');

    %% gft

    gft_adj = U_adj*mem_fn;
    energy_fraction_adj = sqrt(sum(cumsum(abs(gft_adj).^2),2))/norm(gft_adj,'fro');
    r_adj = find(energy_fraction_adj>0.9, 1);
    r_adj = 50;
    
    gft_ha = U_ha*mem_fn;
    energy_fraction_ha = sqrt(sum(cumsum(abs(gft_ha).^2),2))/norm(gft_ha,'fro');
    r_ha = find(energy_fraction_ha>0.9, 1);
    r_ha = 50;
    
    gft_rw = U_rw*mem_fn;
    energy_fraction_rw = sqrt(sum(cumsum(abs(gft_rw).^2),2))/norm(gft_rw,'fro');
    r_rw = find(energy_fraction_rw>0.9, 1);
    r_rw = 50;
    
%     createfigure_gft_bases(1:1:N, [energy_fraction_adj energy_fraction_ha energy_fraction_rw]); 

    %% iterate over signals

    num_trials = 50; % average results for randomized methods

    samples_adj = cell(num_points,1); 
    samples_ha = cell(num_points,1); 
    samples_rw = cell(num_points,1);

    prev_sample_size = 0;

    prev_samples_adj = false(N,1);
    prev_samples_ha = false(N,1);
    prev_samples_rw = false(N,1);

    for i = 1:num_points

        fprintf('******************\n\ni = %d\n\n', i);

        r = sample_size(i);

        % number nodes to be sampled in the current iteration
        sample_size_inc = sample_size(i) - prev_sample_size;

        % sample and reconstruct

        %%%% adj
        [samples_adj{i}, omega_adj,~] = compute_S_L_k_dir(Ln_adj_k_sym, k, sample_size_inc, prev_samples_adj);
        fprintf('number of nodes sampled (proposed) = %d \n\n', sum(samples_adj{i}));
    %     error_adj(i) = reconstruction_pocs_dir(f, fn, samples_adj{i}, Ln_adj, omega_adj);
        error_adj(i,set) = reconstruction_exact_BL_dir(f, fn, samples_adj{i}, V_adj(:,1:r_adj));
        fprintf('classification err adj = %f \n\n', error_adj(i,set));

        %%%% ha
        [samples_ha{i}, omega_ha,~] = compute_S_L_k_dir(Ln_ha_k, k, sample_size_inc, prev_samples_ha);
        fprintf('number of nodes sampled (proposed) = %d \n\n', sum(samples_ha{i}));
    %     error_ha(i) = reconstruction_pocs_dir(f, fn, samples_ha{i}, Ln_ha, omega_ha);
        error_ha(i,set) = reconstruction_exact_BL_dir(f, fn, samples_ha{i}, V_ha(:,1:r_ha));
        fprintf('classification err ha = %f \n\n', error_ha(i,set));

        %%%% rw
        [samples_rw{i}, omega_rw,~] = compute_S_L_k_dir(Ln_rw_k, k, sample_size_inc, prev_samples_rw);
        fprintf('number of nodes sampled (proposed) = %d \n\n', sum(samples_rw{i}));
    %     error_rw(i) = reconstruction_pocs_dir(f, fn, samples_rw{i}, Ln_rw, omega_rw);
        error_rw(i,set) = reconstruction_exact_BL_dir(f, fn, samples_rw{i}, V_rw(:,1:r_rw));
        fprintf('classification err rw = %f \n\n', error_rw(i,set));


        % update variable for next iteration
        prev_sample_size = sample_size(i);

        prev_samples_adj = samples_adj{i};
        prev_samples_ha = samples_ha{i};
        prev_samples_rw = samples_rw{i};
    end


    save('results/usps_results_compare_bases.mat', 'sample_size', 'error_adj' ,'error_ha','error_rw',...
        'N','energy_fraction_adj', 'energy_fraction_ha', 'energy_fraction_rw');


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

