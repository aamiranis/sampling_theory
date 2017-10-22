clear all;

addpath(genpath('util'));

% Define graph size
N = 1000;
K = 8;

A_row = zeros(1,N);

theta = 2*pi/N;
sigma = 2*sin(theta/2);
A_row(1:(K/2+1)) = exp(- (2*sin((0:(K/2))*theta/2)).^2 / sigma^2 );
A_row((N-K/2+1):N) = exp(- (2*sin(((-K/2):-1)*theta/2)).^2 / sigma^2 );

% Define adjacency matrix using random uniformly generated weights
A = zeros(N);
for i=1:N
    A(i,:) = circshift(A_row,[0 i-1]);
end
A = A - diag(diag(A));

% laplacian
[L, Ln, deg] = compute_laplacians(A);

% GFT
[U, lambda] = eig(Ln);
[lambda, idx] = sort(real(diag(lambda)),'ascend');
U = U(:,idx);

save('regular_knn_graph.mat','A', 'deg', 'Ln', 'U', 'lambda');

