function error = rand_sample_reconstruct_uni_dir(f, fn, VR, sample_size)


N = length(f);
S = false(N,1);

S(randsample(N,sample_size)) = 1;


% reconstruct
[error, ~] = reconstruction_exact_BL_dir(f, fn, S, VR);

% f_hat = N * UR * UR(S,:)' * fn(S) / sum(S);
% error = norm(f(~S)-f_hat(~S), 'fro');