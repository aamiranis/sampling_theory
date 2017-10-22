function [ x, lambda, success_flag, log ]...
    = eig_lopcg_cg( A, x, lambda_tol, max_iter, start_precon_tol, precon_tol, max_iter_precon)

%EIG_LOPCG_CG computes the smallest eigenpair of a large sparse symmetric 
% positive definite matrix using conjugate gradients based Rayleigh
% quotient minimization with inbuilt preconditioning based on conjugate
% gradient method again.
% 
% Usage:
% 

    % create matrix vector product handle if it is not
    if isa(A, 'function_handle') 
        A_times_vec = A;
    else
        A_times_vec = @(x) (A*x);
    end

    % calculate rayleigh quotient
    x = x/(sqrt(x'*x));
    v = A_times_vec(x);
    rho = x'*v;
    rho_old = 2*rho;
    x_old = zeros(size(x));
    
    % initializations
    k = 0;
    log = [];
    success_flag = 0;
    
%     while ( rho_old - rho > lambda_tol && k <= max_iter)
    while ( norm(x_old - x) > lambda_tol && k <= max_iter)
       
        k = k+1;
        
        % compute gradient
        g = v - rho*x;
        g_norm = norm(g);
%         g = g / g_norm;
        
        % % % % % % % % % % % % % % % %
        % preconditioning using cg for solving system of equations A*w = g
        % such that w approximates A^(-1) * g;
%         if ( (rho_old - rho < start_precon_tol) && ...
%                 exist('start_precon_tol','var') && ... 
%                 exist('precon_tol','var') && ...
%                 exist('max_iter_precon','var')  )
        if ( (norm(x_old - x) < start_precon_tol) && ...
                        exist('start_precon_tol','var') && ... 
                        exist('precon_tol','var') && ...
                        exist('max_iter_precon','var')  )
            k_cg = 0;
            
            r_cg = g - A_times_vec(g);
            p_cg = r_cg;
            rr_cg_old = r_cg'*r_cg;
            
            while ( (sqrt(rr_cg_old) > precon_tol) && (k_cg <= max_iter_precon) )
                k_cg = k_cg + 1;
                
                A_times_p_int = A_times_vec(p_cg);
                
                alpha_cg = rr_cg_old / (p_cg' * A_times_p_int);
                g = g + alpha_cg * p_cg;
                r_cg = r_cg - alpha_cg * A_times_p_int;
                
                rr_cg_new = r_cg'*r_cg;
                p_cg = r_cg + (rr_cg_new/rr_cg_old) * p_cg;
                
                rr_cg_old = rr_cg_new;
            end
            
        end
        
        g_norm = norm(g);
        g = g / g_norm;
        
        % % % % % % % % % % % % % % % %
        
        % first search direction
        if (k == 1)
            p = zeros(size(x));
        else
            p = p/norm(p);
        end
        
        % store previous rho and x
        rho_old = rho;
        x_old = x;
        
        % Rayleigh-Ritz procedure
        rmat1 = [x -g p]'*[v A_times_vec([-g p])];
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
        
        x = x/(sqrt(x'*x));
        v = A_times_vec(x);
        
        % eigenvalue
        lambda = rho;
        
        % output log
        if(nargout > 3), log = [log; [k, rho, g_norm, rho_old-rho]]; end
    end

    % if iterations converged
    if (k <= max_iter), success_flag = 1; end
    
end
