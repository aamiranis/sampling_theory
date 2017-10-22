function [Lun, Ln, deg] = compute_laplacians(A)

N = length(A);

deg = sum(A,2);
% D = diag(deg);
d = deg;
d = d.^(-1/2);
d(deg==0) = 0;

Dinv = spdiags(d,0,N,N);
Lun = spdiags(deg,0,N,N) - A;
Ln = speye(N) - Dinv*A*Dinv;

% symmetrization of Laplacian
Ln = 0.5*(Ln+Ln');