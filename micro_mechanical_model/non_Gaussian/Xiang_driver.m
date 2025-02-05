clc; clear; close all;
addpath("src/")

% loading data of UT
lambda = importdata("../Hossain_2020/quasistatic_UT/stretch.txt");
P_exp = importdata("../Hossain_2020/quasistatic_UT/stress.txt");
F_list = zeros(3,3,length(lambda));
F_list(1,1,:) = lambda;
F_list(2,2,:) = lambda.^(-0.5);
F_list(3,3,:) = lambda.^(-0.5);


paras0 = [0.01, 100.0];
lb = [0.0, 0.0];
ub = [Inf, Inf];

objectiveFunction = @(paras) objective(paras, F_list, P_exp);

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 40000,  ...
    'TolFun', 1.0e-10, ...
    'TolX', 1.0e-10, ...
    'OptimalityTolerance', 1.0e-10, ...
    'StepTolerance', 1.0e-10, ...
    'Display', 'iter');

[paras, ~] = lsqnonlin( objectiveFunction, paras0, lb, ub, options);
plot_results_Xiang(paras, F_list, P_exp);