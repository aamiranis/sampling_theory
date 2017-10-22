function [error, f_recon] = reconstruction_pocs_dir(f, fn, S_opt, Ln, cutoff)


num_iter = 100;
norm_val = zeros(num_iter,1); % used for checking convergence

% reconstruction using POCS

% approximate low pass filter using SGWT toolbox
filterlen = 10;
alpha = 8;
freq_range = [0 2];
g = @(x)(1./(1+exp(alpha*(x-cutoff))));
c = sgwt_cheby_coeff(g,filterlen,filterlen+1,freq_range);


% initialization
f_init = fn; % noisy samples
f_init(~S_opt,:) = 0;
f_recon = sgwt_cheby_op(f_init,Ln,c,freq_range);

for iter = 1:num_iter % takes fewer iterations
    % projection on C1
    err_s = (f_init-f_recon); 
    err_s(~S_opt,:) = 0; % add the correction only to the known set
    
    % projection on C2
    f_temp = sgwt_cheby_op(f_recon + err_s,Ln,c,freq_range); % err on S approx LP
    
    norm_val(iter) = norm(f_temp-f_recon); % to check convergence
    if (iter > 1 && norm_val(iter) > norm_val(iter-1) ), break; end % avoid divergence
    f_recon = f_temp;
end

f_recon(S_opt,:) = fn(S_opt,:); % observed samples remain unchanged

[~,true_labels] = max(f,[],2);
[~,pred_labels] = max(f_recon,[],2);
error = sum(true_labels(~S_opt) ~= pred_labels(~S_opt))/sum(~S_opt);
% error = norm(f(~S_opt)-f_recon(~S_opt), 'fro')^2;
