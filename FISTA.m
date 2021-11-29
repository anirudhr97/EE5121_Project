function [x,func_eval,numMV] = FISTA(A,b,tau,Gamma,max_iter,type)
x = zeros(size(A,2),1);
numMV = max_iter;
mu = 0;
z = x;
L = norm(A,'fro')^2;
if strcmp(type,'quad_l1')
    for i = 1 : max_iter
        mu = 0.5*(1 + sqrt(1 + 4*mu^2)); gamma = (1-mu)/(mu+1);
        zprev = z;
        z = soft_thresholding((x - (A*x -b)/L), tau/L);
        x = (1-gamma)*z + gamma*zprev;
        func_eval(i) = 0.5*x'*A*x - b'*x + tau*norm(x,1);
    end
elseif strcmp(type,'ls_l1')
    for i = 1 : max_iter
        mu = 0.5*(1 + sqrt(1 + 4*mu^2)); gamma = (1-mu)/(mu+1);
        zprev = z;
        z = soft_thresholding((x - (A'*(A*x -b))/L), tau/L);
        x = (1-gamma)*z + gamma*zprev;
        func_eval(i) = norm(A*x - b,2)^2  + tau*norm(x,1) ;
    end
elseif strcmp(type,'reg_quad_l1')
    for i = 1 : max_iter
        mu = 0.5*(1 + sqrt(1 + 4*mu^2)); gamma = (1-mu)/(mu+1);
        zprev = z;
        z = soft_thresholding((x - ((A+2*gamma*eye(size(A,1)))*x -b)/L), tau/L);
        x = (1-gamma)*z + gamma*zprev;
        func_eval(i) = 0.5*x'*A*x - b'*x + Gamma*norm(x,2) + tau*norm(x,1) ;
    end
end
    
end