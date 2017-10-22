clear all;

addpath(genpath('util'));

% Define graph size
N = 1000;
K = 8;
beta = 0.1;

A_row = zeros(1,N);
A_row(1:(K+1)) = 1;

% Define adjacency matrix
A = zeros(N);
for i=1:N
    A(i,:) = circshift(A_row,[0 i-1-K/2]);
end

% Dont remove self-loops yet, it would help avoid self-loops in random wiring 

% nbr = neighbor
for k = 1:K/2
    % Considering k-th neighbour in clockwise direction
    for i = 1:N-k
        % For row i, look at k-th neighbor, rewire with probability beta
        if(rand(1)<=beta)
            % Look at choices available (self-loop term = 1 avoids rewiring
            % here, later self-loop terms removed).
            choices_available = find(~(A(i,:)>0));
            chosen_node = randsample(choices_available,1);            
            % Perform rewiring
            A(i,i+k) = 0; 
            A(i+k,i) = 0;
            A(i,chosen_node) = 1;
            A(chosen_node,i) = 1;
        end           
    end
end
A = A - diag(diag(A));

% laplacian
[L, Ln, deg] = compute_laplacians(A);

save('small_world_graph_1k.mat','A', 'deg', 'Ln');