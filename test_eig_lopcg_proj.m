function test_eig_lopcg_proj(Ln, S)

k = 8;
N = length(Ln);


% Ln = Ln + 0.1 * speye(length(Ln));

Ln_k = Ln^k;
% S = rand(length(Ln),1) > 0.05;
% S = true(length(Ln),1); S(10) = false;

[y1,s1] = eigs(Ln_k(S,S), 1, 'sm');
s1

% function x = operatorA(x)
%     for i = 1:k
%         x = Ln * x;
%     end
% end
% 
% C = ichol(Ln); prec_fun = @(x)C\(C'\x);
% function x = operatorT(x)
%     for i = 1:k
%         x = prec_fun(x);
%     end
% end

% initial_x = ones(sum(S),1);
% initial_x = initial_x / norm(initial_x);
% [ y2, s2, failure_flag, log ] = eig_lopcg_proj( initial_x, S, @(x)operatorA(x), @(x)operatorT(x), 1e-4, N);
% s2
% plot(log10(log(:,end)));

% initial_x = ones(length(Ln),1);
% initial_x = initial_x / norm(initial_x);
% [ y2, s2 ] = lobpcg( initial_x, @(x)operatorA(x), [], @(x)operatorT(x), 1e-4, 500 );
% y2 = y2(S);
% s2

function result = operatorA(x)
    N = size(Ln,1);
    M = size(x,2);
    result = zeros(N,M);
    
    if (M == 1)
        result(S) = x;
    else
        result(S,:) = x;
    end
    
    for i = 1:k
%         result = (L * result + result) / 3;
        result = Ln * result;
    end
    
    if (M == 1)
        result = result(S);
    else
        result = result(S,:);
    end
end

% function x = operatorT(x)
%     precon_tol = 1e-4;
%     max_iter_precon = 500;
%     
%     k_cg = 0;
% 
%     r_cg = x - operatorA(x);
%     p_cg = r_cg;
%     rr_cg_old = r_cg'*r_cg;
% 
%     while ( (sqrt(rr_cg_old) > precon_tol) && (k_cg <= max_iter_precon) )
%         k_cg = k_cg + 1;
% 
%         A_times_p_int = operatorA(p_cg);
% 
%         alpha_cg = rr_cg_old / (p_cg' * A_times_p_int);
%         x = x + alpha_cg * p_cg;
%         r_cg = r_cg - alpha_cg * A_times_p_int;
% 
%         rr_cg_new = r_cg'*r_cg;
%         p_cg = r_cg + (rr_cg_new/rr_cg_old) * p_cg;
% 
%         rr_cg_old = rr_cg_new;
%     end
% end


% C = ichol(Ln); prec_fun = @(x)C\(C'\x);
% function result = operatorT(x)
%     M = size(x,2);
%     result = zeros(N,M);
%     
%     if (M == 1)
%         result(S) = x;
%     else
%         result(S,:) = x;
%     end
%     
%     for i = 1:k
%         result = prec_fun(result);
%     end
%     
%     if (M == 1)
%         result = result(S);
%     else
%         result = result(S,:);
%     end
% end

initial_x = ones(sum(S),1);
for j = 1:size(initial_x,2)
    initial_x(:,j) = initial_x(:,j)/norm(initial_x(:,j));
end
% [y2, s2, failure_flag, lambda_history, residual_norm_history] = lobpcg(initial_x, @(x)operatorA(x), 1e-4, N);
% plot(log10(residual_norm_history'));
[ y2, s2, failure_flag, log ] = eig_lopcg( initial_x, @(x)operatorA(x), 1e-4, N);
plot(log10(log(:,end)));
s2(1)

figure, plot([y1 y2]);
% plot(abs(y1)-abs(y2));

end