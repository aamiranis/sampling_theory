function error = rand_sample_reconstruct_U_VR_dir(f, fn, U_RV, VR, sample_size)

N = length(f);
num_classes = size(f,2);

% compute the norm of each row of UR
w = sqrt(sum(abs(U_RV).^2, 1));
pdf = w'/sum(w); % pdf over nodes


% plot(pdf)
% pause

% cdf = cumsum(pdf);

S = false(N,1);

S(randsample(N,sample_size,true,pdf)) = 1;

% reconstruct
f_hat = VR * ( U_RV(:,S) * (fn(S,:) ./ repmat(pdf(S),1,num_classes)) ) / sum(S);

[~,true_labels] = max(real(f),[],2);
[~,pred_labels] = max(real(f_hat),[],2);
error = sum(true_labels(~S) ~= pred_labels(~S))/sum(~S);

% error = norm(f(~S)-f_hat(~S), 'fro')^2;
