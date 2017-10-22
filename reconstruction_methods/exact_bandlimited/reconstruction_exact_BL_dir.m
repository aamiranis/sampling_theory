function [error, f_recon] = reconstruction_exact_BL_dir(f, fn, S, VR)

V_SR = VR(S,:);

f_recon = VR * (V_SR \ fn(S,:));
% f_recon = VR * (pinv(V_SR) * fn(S,:));
f_recon(S,:) = fn(S,:); % observed samples are unchanged

[~,true_labels] = max(real(f),[],2);
[~,pred_labels] = max(real(f_recon),[],2);
error = sum(true_labels(~S) ~= pred_labels(~S))/sum(~S);
% error = norm(f(~S)-f_recon(~S),'fro')^2;