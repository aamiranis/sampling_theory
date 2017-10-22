clear all;

addpath(genpath('../util'));

n = 100000;
% p = 10/n;
p = 0.01;

rng('default');

G = sprand(n,n,p);
G = G > 0;
G = triu(G,1);
G = G + G';
Sum = sum(G,2);
ind = find(Sum ~=0);
G = G(ind,ind);

A = G;

% A = sparse(A);

% laplacian
[L, Ln, deg] = compute_laplacians(A);

% save('erdos_renyi_graph_50k.mat','A', 'deg', 'Ln');