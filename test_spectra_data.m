%% Testing the working of iiCG, FISTA and ISTA on the Spectra Problem
% Closing all open figures
close all;

%{
-- If you do not have the spectra dataset on your MATLAB version, use the 
'spectra.mat' file given in this repository. Use the following command 
to import the dataset into your MATLAB workspace:
$ spectra = importdata('spectra.mat')

-- This code can be executed only after the dataset has been imported 
into your workspace.

-- Use the following command to get an idea on what the dataset is about:
$ spectra.Description
%}

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
max_iter = length(fvals_iicg)+15;
[~,func_eval_fista,~,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Running ISTA
tau = tau(1);
max_iter = length(fvals_iicg)+15;
[~,func_eval_ista,~,~]  =  ISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Finding the lowest function value among the 3 methods.
Fstar = min([min(fvals_iicg), min(func_eval_ista), min(func_eval_fista)]);

% Calculating tolerances
tol_iiCG = (fvals_iicg - Fstar)/abs(Fstar);
tol_fista = (func_eval_fista - Fstar)/abs(Fstar);
tol_ista = (func_eval_ista - Fstar)/abs(Fstar);

% Plotting the tolerance for each of ISTA, FISTA and iiCG
figure;
plot(tol_ista, 'r', 'LineWidth', 1.5);
hold on;
plot(tol_fista, 'k', 'LineWidth', 1.5);
hold on;
plot(tol_iiCG, 'b', 'LineWidth', 1.5);
legend('ISTA', 'FISTA', 'iiCG');
grid on;
xlabel('Number of Matrix-Vector Products', 'Interpreter','latex', 'FontSize', 13)
ylabel('Tolerance ($\frac{F(x^t)-F^*}{|F^*|}$)', 'Interpreter','latex', 'FontSize', 13)
title('Plot of Tolerance vs Matrix-Vector Products for Spectra Problem')
print('Plots/tolvsMV','-dpng');


% Plotting with logarithmic y axis
figure;
semilogy(tol_ista, 'r', 'LineWidth', 1.5);
hold on;
semilogy(tol_fista, 'k', 'LineWidth', 1.5);
hold on;
semilogy(tol_iiCG, 'b', 'LineWidth', 1.5);
legend('ISTA', 'FISTA', 'iiCG');
grid on;
ylim([1e-12 inf]);
xlabel('Number of Matrix-Vector Products', 'Interpreter','latex', 'FontSize', 13)
ylabel('Tolerance ($\frac{F(x^t)-F^*}{|F^*|}$)', 'Interpreter','latex', 'FontSize', 13)
title('Plot of Tolerance vs Matrix-Vector Products for Spectra Problem')
print('Plots/tolvsMV_log','-dpng');


% Plotting CG Move Count for iiCG
figure;
plot(out2.CGmvcount, 'LineWidth', 1.5);
grid on;
xlabel('Iterations')
ylabel('CG Move Count')
title('Plot of Number of CG moves vs Iterations for Spectra Problem')
print('Plots/CG','-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting for more number of iterations

% Running FISTA
tau = tau(1);
max_iter = 800;
[x_fista,func_eval_fista,~,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Running ISTA
tau = tau(1);
max_iter = 800;
[x_ista,func_eval_ista,~,~]  =  ISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Finding the lowest function value among the 3 methods.
Fstar = min([min(fvals_iicg), min(func_eval_ista), min(func_eval_fista)]);

% Calculating tolerances
tol_iiCG = (fvals_iicg - Fstar)/abs(Fstar);
tol_fista = (func_eval_fista - Fstar)/abs(Fstar);
tol_ista = (func_eval_ista - Fstar)/abs(Fstar);

% Plotting the tolerance for each of ISTA, FISTA and iiCG
figure;
plot(tol_ista, 'r', 'LineWidth', 1.5);
hold on;
plot(tol_fista, 'k', 'LineWidth', 1.5);
hold on;
plot(tol_iiCG, 'b', 'LineWidth', 1.5);
legend('ISTA', 'FISTA', 'iiCG');
grid on;
xlabel('Number of Matrix-Vector Products', 'Interpreter','latex', 'FontSize', 13)
ylabel('Tolerance ($\frac{F(x^t)-F^*}{|F^*|}$)', 'Interpreter','latex', 'FontSize', 13)
title('Plot of Tolerance vs Matrix-Vector Products for Spectra Problem')
print('Plots/tolvsMVlll','-dpng');


% Plotting with logarithmic y axis
figure;
semilogy(tol_ista, 'r', 'LineWidth', 1.5);
hold on;
semilogy(tol_fista, 'k', 'LineWidth', 1.5);
hold on;
semilogy(tol_iiCG, 'b', 'LineWidth', 1.5);
legend('ISTA', 'FISTA', 'iiCG');
grid on;
ylim([1e-12 inf]);
xlabel('Number of Matrix-Vector Products', 'Interpreter','latex', 'FontSize', 13)
ylabel('Tolerance ($\frac{F(x^t)-F^*}{|F^*|}$)', 'Interpreter','latex', 'FontSize', 13)
title('Plot of Tolerance vs Matrix-Vector Products for Spectra Problem')
print('Plots/tolvsMV_loglll','-dpng');
