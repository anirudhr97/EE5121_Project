function [x,func_eval,numMV,xHist,zHist] = FISTA(A,b,tau,Gamma,max_iter,type,opts)
%{
Function to implement the Fast Iterative Shrinkage Thresholding Algorithm
for 3 cases of objective function:
                1)Quadratic       + L1 regularization
                2)Least Squares   + L1 regularization
                3)Quadratic       + L1 regularization     + L2 regularization
Though these can be reduced to the same form, all three forms are allowed
for convenience
'Quadratic' here refers to the form x^TAx - bx
%}
arguments
   A double
   b double
   tau double % coefficient for l1 regulariser
   Gamma double % applies only to (3), coefficient for L2 regulariser
   max_iter {mustBeInteger}
   type {mustBeMember(type,{'quad_l1','ls_l1','reg_quad_l1'})}
   opts.x0 double = zeros(size(A,2),1)
end
x = opts.x0;
numMV = max_iter;
func_eval = zeros(1,max_iter+1);
xHist = zeros(max_iter+1, length(x));
zHist = zeros(max_iter+1, length(x));
xHist(1,:) = x;
zHist(1,:) = x;
mu = 0; % initial momentum acceleration parameter
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