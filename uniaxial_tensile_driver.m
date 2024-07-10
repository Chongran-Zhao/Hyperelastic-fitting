clc; clear; close all;
addpath("src/")
% loading data
lambda = importdata("./Treloar-UT/strain.txt");
P_exp = importdata("./Treloar-UT/stress.txt");

Ft = zeros(3,3,length(lambda));
Ft(1,1,:) = lambda;
Ft(2,2,:) = lambda.^(-0.5);
Ft(3,3,:) = lambda.^(-0.5);

mu = [1.0, 1.0, 1.0];
alpha = [1.0, 1.0, 1.0];

[paras0, lb, ub, num] = array_to_paras(mu, alpha);
objectiveFunction = @(paras) objective(paras, Ft, P_exp, num);
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'interior-point', ...
    'MaxIterations', 1000, ...
    'Display', 'iter', ...
    'TolFun', 1e-16, ...
    'TolX', 1e-16);
[paras, ~] = lsqnonlin( objectiveFunction, paras0, lb, ub, options);

plot_result(paras, num, Ft, P_exp);
print(gcf, '-djpeg', 'fig.jpg');