function [x,func_eval,numMV,xHist,zHist] = FISTA(A,b,tau,Gamma,max_iter,type)
%{
Function to implement the Fast Iterative Shrinkage Thresholding Algorithm
for 3 cases:    Quadratic       + L1 regularization
                Least Squares   + L1 regularization
                Quadratic       + L1 regularization     + L2 regularization
%}
x = zeros(size(A,2),1);
numMV = max_iter;
func_eval = zeros(max_iter+1, 1);
xHist = zeros(max_iter+1, length(x));
zHist = zeros(max_iter+1, length(x));
xHist(1,:) = x;
zHist(1,:) = x;
mu = 0;
z = x;
L = norm(A,'fro')^2;
if strcmp(type,'quad_l1')
    func_eval(1) = 0.5*x'*A*x - b'*x + tau*norm(x,1);
    for i = 1 : max_iter
        mu = 0.5*(1 + sqrt(1 + 4*mu^2)); gamma = (1-mu)/(mu+1);
        zprev = z;
        z = soft_thresholding((x - (A*x -b)/L), tau/L);
        x = (1-gamma)*z + gamma*zprev;
        func_eval(i+1) = 0.5*x'*A*x - b'*x + tau*norm(x,1);
        xHist(i+1,:) = x;
        zHist(i+1,:) = z;
    end
elseif strcmp(type,'ls_l1')
    func_eval(1) = norm(A*x - b,2)^2  + tau*norm(x,1);
    for i = 1 : max_iter
        mu = 0.5*(1 + sqrt(1 + 4*mu^2)); gamma = (1-mu)/(mu+1);
        zprev = z;
        z = soft_thresholding((x - (A'*(A*x -b))/L), tau/L);
        x = (1-gamma)*z + gamma*zprev;
        func_eval(i+1) = norm(A*x - b,2)^2  + tau*norm(x,1);
        xHist(i+1,:) = x;
        zHist(i+1,:) = z;
    end
elseif strcmp(type,'reg_quad_l1')
    func_eval(1) = 0.5*x'*A*x - b'*x + Gamma*norm(x,2) + tau*norm(x,1);
    for i = 1 : max_iter
        mu = 0.5*(1 + sqrt(1 + 4*mu^2)); gamma = (1-mu)/(mu+1);
        zprev = z;
        z = soft_thresholding((x - ((A+2*gamma*eye(size(A,1)))*x -b)/L), tau/L);
        x = (1-gamma)*z + gamma*zprev;
        func_eval(i+1) = 0.5*x'*A*x - b'*x + Gamma*norm(x,2) + tau*norm(x,1);
        xHist(i+1,:) = x;
        zHist(i+1,:) = z;
    end
end
end