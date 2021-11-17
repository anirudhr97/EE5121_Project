% Create 2D/3D visualisations of steps in the algorithm, and in FISTA
clear;
close all;

% rng(3); % seed for reproducibility

% 2D problem
pdim = 2;
M = rand(pdim,pdim);
A  = M * M'; % sym pos def matrix
% b  = rand(pdim,1);
b = [0.5; 0.5];
tau = 0.5;

% Conditioning Required?
Ax = @(x) A*x;

% Run algo
problem.Ax = Ax;
problem.b =b;
problem.tau = tau*ones(pdim,1);
opts.x_0 = [3;2];
% [out1] = alg_ql1(problem,opts);
[out1,out2,xPrev] = alg_ql1(problem,opts); % get history of x

disp('Algo Status');
disp(out1.algStatus);


% Plot contours of the function

% ql1 function F (for contour plotting)
F_2D = @ (X) sum(X.*(A*X),1) - b'*X  + tau * vecnorm(X,1); % X is a d x N input (N sample points)
N = 50; % discretization no
viewbound = 1.5*max(max(abs(xPrev),[],'all'),1);
x1 = linspace(-viewbound,viewbound,N);
x2 = linspace(-viewbound,viewbound,N);
[X1,X2] = meshgrid(x1,x2);
X_wrap = reshape(cat(3,X1,X2),N*N,2);
Z = reshape(F_2D(X_wrap'),N,N);
contourf(X1,X2,Z,10);
grid('on');
hold on;
% Plot history of x, along with psi and phi
plot(xPrev(1,:),xPrev(2,:),'xr-.');
title('Path of solution');

% 3D problem (will be difficult to show contours, but nice to show sparsity, and exploring a plane using CG)
