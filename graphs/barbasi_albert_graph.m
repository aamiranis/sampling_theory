% This file generates a scale-free graph in accordance with the
% Barbasi-Albert model.

clear all;

addpath(genpath('util'));

rng('default');

N = 1000;
m = 4;
d = 4;

G = ones(m);
G = G - diag(diag(G));

A = zeros(N);
A(1:m,1:m) = G;

for node = (m+1) : N
    degrees = sum(A(1:(node-1), 1:(node-1)),2);
    neighbors = randsample(1:(node-1), d, true, degrees);
    A(node, neighbors) = 1;
    A(neighbors, node) = 1;
end

% laplacian
[L, Ln, deg] = compute_laplacians(A);

% GFT
[U, lambda] = eig(Ln);
[lambda, idx] = sort(real(diag(lambda)),'ascend');
U = U(:,idx);

save('barbasi_albert_graph.mat','A', 'deg', 'Ln', 'U', 'lambda');