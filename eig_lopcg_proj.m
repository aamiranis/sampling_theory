function [ x, lambda, success_flag, log ]...
    = eig_lopcg_proj( initial_x, proj_set, operatorA, operatorT, tol, max_iter)

%EIG_LOPCG_CG computes the smallest eigenpair of a large sparse symmetric 
% positive definite matrix using conjugate gradients based Rayleigh
% quotient minimization with inbuilt preconditioning based on conjugate
% gradient method again.
% 
% Usage:
% 
    
    x = zeros(size(proj_set,1),1);
    x(proj_set) = initial_x;

    % calculate rayleigh quotient
    x = x/norm(x);
    v = operatorA(x);
    
    rho = x'*v;
    
    x_old = zeros(size(x));
    
    % initializations
    k = 0;
    log = [];
    success_flag = 0;
    
    while ( norm(x_old - x) > tol && k <= max_iter)
       
        k = k+1;
        
        % compute gradient
        g = v - rho*x;
        g(~proj_set) = 0;
%         g = g / norm(g);
        
        % Precondition g
%         g = operatorT(g);

        g(~proj_set) = 0;
%         g = g / norm(g);
        
        % % % % % % % % % % % % % % % %
        
        % first search direction
        if (k == 1)
            p = zeros(size(x));
        else
%             p = p/norm(p);
        end
        
        % store previous rho and x
        rho_old = rho;
        x_old = x;
        
        % Rayleigh-Ritz procedure
%         [Q,~] = qr([x g p]);
%         x = Q(:,1); g = Q(:,2); p = Q(:,3);
        rmat1 = [x -g p]'*[v operatorA([-g p])];
        rmat2 = [x -g p]'*[x -g p];
        rmat1 = (rmat1+rmat1')/2;
        rmat2 = (rmat2+rmat2')/2;
        [rvec, rval] = eig(rmat1, rmat2);
        [rho, index] = min(real(diag(rval)));
        coeffs = rvec(:,index);
        
        % new search direction
        p = [-g p]*coeffs(2:end);
        
        % new eigenvector
        x = coeffs(1)*x + p;
       
        % Projection step
        x(~proj_set) = 0;        
        x = x/norm(x);
        
        v = operatorA(x);
        
        % eigenvalue
        lambda = rho;
        
        % output log
        if(nargout > 3), log = [log; [k, rho, rho_old-rho, norm(x_old-x)]]; end
    end

    % if iterations converged
    if (k <= max_iter), success_flag = 1; end
    
    x = x(proj_set);
    
end
