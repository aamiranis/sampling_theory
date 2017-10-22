function error = rand_sample_reconstruct_U_VR(f, fn, UR, sample_size)

N = length(f);

% compute the norm of each row of UR
w = sqrt(sum(abs(UR).^2, 2));
pdf = w/sum(w); % pdf over nodes

% plot(pdf)
% pause

% cdf = cumsum(pdf);

S = false(N,1);

% sample with the above pdf
% while (sum(S) < sample_size)
%     y = rand(1);
%     idx = find(cdf<=y, 1, 'last');
%     S(idx) = 1;
% end
S(randsample(N,sample_size,true,pdf)) = 1;

% reconstruct using method presented in paper
% f_hat = UR * ( UR(S,:)' * (fn(S)./w(S)) ) / sum(S);
% 
% error = norm(f(~S)-f_hat(~S), 'fro')^2;

% reconstruct using exact BL
USR = UR(S,:);

f_recon = UR * (USR \ fn(S));
f_recon(S) = fn(S); % observed samples are unchanged

error = norm(f(~S)-f_recon(~S),'fro')^2;
