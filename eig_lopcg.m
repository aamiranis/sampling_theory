function [ x, lambda, success_flag, log ]...
    = eig_lopcg( x, operatorA, tol, max_iter)

    % calculate rayleigh quotient
    x = x/(sqrt(x'*x));
    Ax = operatorA(x);
    rho = x'*Ax;
    
    % initializations
    k = 0;
    log = [];
    success_flag = 0;
    
    g_norm = 1;
    
    while ( g_norm > tol && k <= max_iter)
       
        k = k+1;
        
        % compute gradient
        g = Ax - rho*x;
%         g = Ax;
        g_norm = norm(g);
        g = g/g_norm;
    
        
        % first search direction
        if (k == 1)
            p = zeros(size(x));
        else
            p = p/norm(p);
        end
        
        % Rayleigh-Ritz procedure
%         [Q,~] = qr([x g p]);
%         x = Q(:,1); g = Q(:,2); p = Q(:,3);
        rmat1 = [x -g p]'*[Ax operatorA([-g p])];
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
        x = x/norm(x);
        
        Ax = operatorA(x);
        
        % eigenvalue
        lambda = rho;
        
        % output log
        if(nargout > 3), log = [log; [k, rho, g_norm]]; end
    end

    % if iterations converged
    if (k <= max_iter), success_flag = 1; end
    
end
