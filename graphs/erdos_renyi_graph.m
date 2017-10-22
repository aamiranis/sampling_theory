clear all;

addpath(genpath('util'));

n = 1000;
p = 0.1

rng('default');

G = rand(n);
G = G < p;
G = triu(G,1);
G = G + G';
Sum = sum(G,2);
ind = find(Sum ~=0);
G = G(ind,ind);

A = G;

% laplacian
[L, Ln, deg] = compute_laplacians(A);

% GFT
[U, lambda] = eig(Ln);
[lambda, idx] = sort(real(diag(lambda)),'ascend');
U = U(:,idx);

save('erdos_renyi_graph_10.mat','A', 'deg', 'Ln', 'U', 'lambda');