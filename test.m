% 
rng(33);
Gamma = 0.1;
B = randn(100,200); 
y = 200*randn(100,1);

A = B'*B + 2*Gamma*eye(size(B,2));
b = B'*y;

tau = 1* ones(size(B,2),1);

problem.Ax = @(x) A*x;
problem.tau = tau;
problem.b = b;

[out1,out2]  = alg_ql1(problem);
fvals_iicg = out2.fValues;

% FISTA
tau = tau(1);
disp('FISTA:');

max_iter = 1500;
[x,func_eval,numMV,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');
disp(func_eval(end));

% Fstar = min([fvals_iicg(end),func_eval(end)]);
% tol = 

figure;
plot(func_eval);
figure
plot(fvals_iicg);