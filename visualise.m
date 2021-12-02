% Create 2D visualisations of steps in the algorithm, and in FISTA
clear;
close all;

rng(5);

% Define 2D problem
pdim = 2; % dimension of problem

% Random problem
% M = rand(pdim,pdim);
% A  = M * M'; % sym pos def matrixplot(xHist_fista(:,1),xHist_fista(:,2),'og-.');
% [E,D] = eig(A);
% E = [1 1; 1 -1];
% D = diag([1 3]); % make the quadratic part almost circular
% A = E*D*E';
% Handcrafted choice
A = diag([1.5 0.8]);
b = [1; 0];
tau = 0.5; % regularisation param for L1 regularisation

% Function handle to return matrix-vector product
Ax = @(x) A*x;

% Starting position
% x0=[0;0];
x0 = [1;-0.5];

%% iiCG
problem.Ax = Ax;
problem.b =b;
problem.tau = tau*ones(pdim,1);
opts.x_0 = x0;
opts.separate = true;
% [out1] = alg_ql1(problem,opts);
[out1,out2,xPrev] = alg_ql1(problem,opts); % get history of x
disp('Algo Status');
disp(out1.algStatus);

% These are run for the same number of matrix-vector products
%% FISTA (with same initialisation)
disp('FISTA:');
[x_fista,func_eval_fista,numMV_fista,xHist_fista,zHist_fista] = FISTA(A,b,tau,NaN,out1.numA,'quad_l1','x0',opts.x_0);
disp(func_eval_fista(end));

%% ISTA
% disp('ISTA:');
[x_ista,func_eval_ista,numMV_ista,xHist_ista] = ISTA(A,b,tau,NaN,out1.numA,'quad_l1','x0',opts.x_0);
disp(func_eval_ista(end));

%% Plot contours of the function
center = out1.x; % make sure optimal in the center of the plot

% ql1 function F (for contour plotting)
F_2D = @ (X) 0.5*sum(X.*(A*X),1) - b'*X  + tau * vecnorm(X,1); % X is a d x N input (N sample points)
N = 100; % discretization no
% find bounds for plotting
x1min = min([xPrev(1,:).';xHist_ista(:,1);xHist_fista(:,1)]);
x1max = max([xPrev(1,:).';xHist_ista(:,1);xHist_fista(:,1)]);
x2min = min([xPrev(2,:).';xHist_ista(:,2);xHist_fista(:,2)]);
x2max = max([xPrev(2,:).';xHist_ista(:,2);xHist_fista(:,2)]);
x1min = min(center(1)-1,x1min);x1max = max(center(1)+1,x1max);
x2min = min(center(2)-1,x2min); x2max = max(center(2)+1,x2max);


x1 = linspace(x1min-0.1,x1max+0.1,N);
x2 = linspace(x2min-0.1,x2max+0.1,N);
[X1,X2] = meshgrid(x1,x2);
X_wrap = reshape(cat(3,X1,X2),N*N,2);
Z = reshape(F_2D(X_wrap'),N,N);
contour(X1,X2,Z,10,'DisplayName','');
grid('on');
hold on;
% Plot history of x fpr iiCG
h1 = plot(xPrev(1,:),xPrev(2,:),'or-','DisplayName','iiCG');
% for FISTA
h2 = plot(xHist_fista(:,1),xHist_fista(:,2),'og-.','DisplayName','FISTA');
% plot(xHist_ista(:,1),xHist_ista(:,2),'oc-.');
hold off;
xline(0,'DisplayName','');
yline(0,'DisplayName','');
legend([h1 h2],'Location','northwest');
xlabel('x_1');
ylabel('x_2');
title('Solution Path');