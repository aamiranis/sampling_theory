function [error, f_recon] = reconstruction_exact_BL(f, fn, S, UR)

USR = UR(S,:);

f_recon = UR * (USR \ fn(S));
f_recon(S) = fn(S); % observed samples are unchanged

error = norm(f(~S)-f_recon(~S),'fro')^2;