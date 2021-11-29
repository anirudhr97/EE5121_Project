%% Testing the working of iiCG, FISTA and ISTA on the Spectra Problem
% Closing all open figures
% close all;

% Obtaining B and y values from the spectra dataset
Gamma = 0.1;
B = spectra.NIR;        % 60 x 401
y = spectra.octane;     % 60 x 1

% Forming A and b corresponding to the quadratic L1 problem form of the 
% Spectra Problem.
A = B'*B + 2*Gamma*eye(size(B,2));
b = B'*y;
tau = 1* ones(size(B,2),1);

% Running iiCG
problem.Ax = @(x) A*x;
problem.tau = tau;
problem.b = b;
[out1,out2]  = alg_ql1(problem);
fvals_iicg = out2.fValues;

% Running FISTA
tau = tau(1);
max_iter = 700;
[x_fista,func_eval_fista,~,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Running ISTA
tau = tau(1);
max_iter = 700;
[x_ista,func_eval_ista,~,~]  =  ISTA_final(A,b,tau,Gamma,max_iter,'quad_l1');

% Finding the lowest function value among the 3 methods.
Fstar = min([fvals_iicg(end),func_eval_ista(end),func_eval_fista(end)]);

% Calculating tolerances
tol_iiCG = (fvals_iicg-Fstar)/abs(Fstar);
tol_fista = (func_eval_fista-Fstar)/abs(Fstar);
tol_ista = (func_eval_ista-Fstar)/abs(Fstar);

% Plotting the tolerance for each of ISTA, FISTA and iiCG
figure;
plot(func_eval_fista);
hold on;
plot(func_eval_ista);
hold on;
plot(fvals_iicg);
legend('FISTA', 'ISTA', 'iiCG');
grid on;
xlabel('Number of Matrix-Vector Products')
ylabel('Tolerance')
title('Plot of Tolerance vs Matrix-Vector Products for Spectra Problem')
savefig('tolvsMV.fig');

% Plotting with logarithmic x axis
figure;
semilogx(func_eval_fista);
hold on;
semilogx(func_eval_ista);
hold on;
semilogx(fvals_iicg);
legend('FISTA', 'ISTA', 'iiCG');
grid on;
xlabel('Number of Matrix-Vector Products')
ylabel('Tolerance')
title('Plot of Tolerance vs Matrix-Vector Products for Spectra Problem')
savefig('tolvsMV_logx.fig');
