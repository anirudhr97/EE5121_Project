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
[x_fista,func_eval_fista,~,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');

% Running ISTA
tau = tau(1);
max_iter = length(fvals_iicg)+100;
[x_ista,func_eval_ista,~,~]  =  ISTA_final(A,b,tau,Gamma,max_iter,'quad_l1');

% Finding the lowest function value among the 3 methods.
Fstar = min([min(fvals_iicg), min(func_eval_ista), min(func_eval_fista)]);

% Calculating tolerances
tol_iiCG = (fvals_iicg - Fstar)/abs(Fstar);
tol_fista = (func_eval_fista - Fstar)/abs(Fstar);
tol_ista = (func_eval_ista - Fstar)/abs(Fstar);

% Plotting the tolerance for each of ISTA, FISTA and iiCG
figure;
plot(tol_ista, 'r');
hold on;
plot(tol_fista, 'g');
hold on;
plot(tol_iiCG, 'b');
legend('ISTA', 'FISTA', 'iiCG');
grid on;
xlabel('Number of Matrix-Vector Products', 'Interpreter','latex', 'FontSize', 13)
ylabel('Tolerance ($\frac{F(x^t)-F^*}{|F^*|}$)', 'Interpreter','latex', 'FontSize', 13)
title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
savefig('Plots/tolvsMV_rand.fig');

% Plotting with logarithmic x axis
figure;
semilogy(tol_ista, 'r');
hold on;
semilogy(tol_fista, 'g');
hold on;
semilogy(tol_iiCG, 'b');
legend('ISTA', 'FISTA', 'iiCG');
grid on;
xlabel('Number of Matrix-Vector Products', 'Interpreter','latex', 'FontSize', 13)
ylabel('Tolerance ($\frac{F(x^t)-F^*}{|F^*|}$)', 'Interpreter','latex', 'FontSize', 13)
title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
savefig('Plots/tolvsMV_rand_logy.fig');

% Plotting CG Move Count for iiCG
figure;
plot(out2.CGmvcount);
grid on;
xlabel('Iterations')
ylabel('CG Move Count')
title('Plot of Number of CG moves vs Iterations for Random Matrix Problem')
savefig('Plots/CGvsiter_rand.fig');

%% Commented out code
% % Plotting the tolerance for each of ISTA, FISTA and iiCG
% figure;
% plot(func_eval_ista, 'r');
% hold on;
% plot(func_eval_fista, 'g');
% hold on;
% plot(fvals_iicg, 'b');
% legend('ISTA', 'FISTA', 'iiCG');
% grid on;
% xlabel('Number of Matrix-Vector Products')
% ylabel('Tolerance')
% title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
% savefig('tolvsMV_rand.fig');
% 
% % Plotting with logarithmic x axis
% figure;
% semilogx(func_eval_ista, 'r');
% hold on;
% semilogx(func_eval_fista, 'g');
% hold on;
% semilogx(fvals_iicg, 'b');
% legend('ISTA', 'FISTA', 'iiCG');
% grid on;
% xlabel('Number of Matrix-Vector Products')
% ylabel('Tolerance')
% title('Plot of Tolerance vs Matrix-Vector Products for Random Matrix Problem')
% savefig('tolvsMV_rand_logx.fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 
% rng(33);
% Gamma = 0.1;
% B = randn(100,200); 
% y = 200*randn(100,1);
% 
% A = B'*B + 2*Gamma*eye(size(B,2));
% b = B'*y;
% 
% tau = 1* ones(size(B,2),1);
% 
% problem.Ax = @(x) A*x;
% problem.tau = tau;
% problem.b = b;
% 
% [out1,out2]  = alg_ql1(problem);
% fvals_iicg = out2.fValues;
% 
% % FISTA
% tau = tau(1);
% disp('FISTA:');
% 
% max_iter = 766;
% [x,func_eval,numMV,~,~]  =  FISTA(A,b,tau,Gamma,max_iter,'quad_l1');
% disp(func_eval(end));
% 
% % Fstar = min([fvals_iicg(end),func_eval(end)]);
% % tol = 
% close all
% figure;
% plot(func_eval);
% hold on;
% plot(fvals_iicg);
% 
% % no of CG cycless
% plot(1:length(out2.CGmvcount),out2.CGmvcount)
