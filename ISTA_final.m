function [x,func_eval,numMV,xHist] = ISTA_final(A,b,tau,gamma,max_iter,type,opts)
%{
Function to implement the Iterative Shrinkage Thresholding Algorithm
for 3 cases:    Quadratic       + L1 regularization
                Least Squares   + L1 regularization
                Quadratic       + L1 regularization     + L2 regularization
%}
arguments
   A double
   b double
   tau double
   gamma double
   max_iter {mustBeInteger}
   type {mustBeMember(type,{'quad_l1','ls_l1','reg_quad_l1'})}
   opts.x0 double = zeros(size(A,2),1)
end
x = opts.x0;
numMV = max_iter;
func_eval = zeros(1,max_iter+1);
xHist = zeros(max_iter+1, length(x));
xHist(1,:) = x;
L = norm(A,'fro')^2;
if strcmp(type,'quad_l1')
    func_eval(1) = 0.5*x'*A*x - b'*x + tau*norm(x,1);
    for i = 1 : max_iter
        x = soft_thresholding((x - (A*x -b)/L), tau/L);
        xHist(i+1,:) = x;
        func_eval(i+1) = 0.5*x'*A*x - b'*x + tau*norm(x,1);
    end
elseif strcmp(type,'reg_quad_l1')
    func_eval(1) = 0.5*x'*A*x - b'*x + gamma*norm(x,2) + tau*norm(x,1) ;
    for i = 1 : max_iter
        x = soft_thresholding((x - ((A+2*gamma*eye(size(A,1)))*x -b)/L), tau/L);
        xHist(i+1,:) = x;
        func_eval(i+1) = 0.5*x'*A*x - b'*x + gamma*norm(x,2) + tau*norm(x,1) ;
    end
elseif strcmp(type,'ls_l1')
    func_eval(1) = norm(A*x - b,2)^2  + tau*norm(x,1) ;
    for i = 1 : max_iter
        x = soft_thresholding((x - (A'*(A*x -b))/L), tau/L);
        xHist(i+1,:) = x;
        func_eval(i+1) = norm(A*x - b,2)^2  + tau*norm(x,1) ;
    end
end

end




