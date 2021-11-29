
function [x,func_eval] = ISTA_final(A,b,tau,gamma,max_iter,type)
x = zeros(size(A,2),1); 
L = norm(A,'fro')^2;
if strcmp(type,'quad_l1')
    for i = 1 : max_iter
        x = soft_thresholding((x - (A*x -b)/L), tau/L);
        func_eval(i) = 0.5*x'*A*x - b'*x + tau*norm(x,1);
    end
elseif strcmp(type,'reg_quad_l1')
    for i = 1 : max_iter
        x = soft_thresholding((x - ((A+2*gamma*eye(size(A,1)))*x -b)/L), tau/L);
        func_eval(i) = 0.5*x'*A*x - b'*x + gamma*norm(x,2) + tau*norm(x,1) ;
    end
elseif strcmp(type,'ls_l1')
    for i = 1 : max_iter
        x = soft_thresholding((x - (A'*(A*x -b))/L), tau/L);
        func_eval(i) = norm(A*x - b,2)^2  + tau*norm(x,1) ;
    end
end

end




