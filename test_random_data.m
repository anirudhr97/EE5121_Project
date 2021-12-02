%% Testing the working of iiCG, FISTA and ISTA on the Random Matrix Problem
% Closing all open figures
close all;

% Constructing B and y as random matrices from the normal distribution
rng(3);            % For repeatability
Gamma = 0.1;
B = randn(100,200);
y = 200*randn(100,1);

% Forming A and b corresponding to the quadratic L1 problem form of the 
% Random Matrix Problem.
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
max_iter = length(fvals_iicg)+100;
[~,func_eval_fista,~,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Running ISTA
tau = tau(1);
max_iter = length(fvals_iicg)+100;
[~,func_eval_ista,~,~]  =  ISTA_final(A,b,tau,Gamma,max_iter,'quad_l1');

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
title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
print('Plots/tolvsMV_rand','-dpng');

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
title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
print('Plots/tolvsMV_log_rand','-dpng');

% Problem specific calculations to plot the no. of CG moves as a bar graph
% First few values are omitted to give a better number to plot a bar graph
% They were anyways constant. So, nothing interesting...
arrr = out2.CGmvcount(5:end);
p = length(arrr);
a = sum(reshape(arrr,12,p/12));
labels = zeros(1,22);
for i=1:22
    labels(i) = strcat(string(12*i-11), '-', string(12*i));
end

% Plotting CG Move Count for iiCG
figure;
bar(a, 'LineWidth', 1.5);
set(gca, 'XTick', 1:22, 'XTickLabel',labels);
grid on;
xlabel('Iterations')
ylabel('CG Move Count')
title('Plot of Number of CG moves vs Iterations for Random Matrix Problem')
print('Plots/CG_rand','-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting for more number of iterations

% Running FISTA
tau = tau(1);
max_iter = 1500;
[x_fista,func_eval_fista,~,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Running ISTA
tau = tau(1);
max_iter = 1500;
[x_ista,func_eval_ista,~,~]  =  ISTA_final(A,b,tau,Gamma,max_iter,'quad_l1');

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
title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
print('Plots/tolvsMV_randlll','-dpng');

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
title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
print('Plots/tolvsMV_log_randlll','-dpng');
