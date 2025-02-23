clc; clear; close all;
addpath("src/")

% loading data of UT
% data_1 = importdata("../Treloar_1944/UT/stress_stretch.txt");
% data_1 = importdata("../Kawabata_1981/UT/stress_stretch.txt");
data_1 = importdata("../Meunier_2008/UT/stress_stretch.txt");
P_exp_1 = data_1(:,1);
lambda_1 = data_1(:,2);
F_list_1 = zeros(3,3,length(lambda_1));
F_list_1(1,1,:) = lambda_1;
F_list_1(2,2,:) = lambda_1.^(-0.5);
F_list_1(3,3,:) = lambda_1.^(-0.5);

% loading data of ET
% data_2 = importdata("../Treloar_1944/ET/stress_stretch.txt");
% data_2 = importdata("../Kawabata_1981/ET/stress_stretch.txt");
data_2 = importdata("../Meunier_2008/ET/stress_stretch.txt");
P_exp_2 = data_2(:,1);
lambda_2 = data_2(:,2);
F_list_2 = zeros(3,3,length(lambda_2));
F_list_2(1,1,:) = lambda_2;
F_list_2(2,2,:) = lambda_2;
F_list_2(3,3,:) = lambda_2.^(-2.0);

% loading data of PS
% data_3 = importdata("../Treloar_1944/PS/stress_stretch.txt");
% data_3 = importdata("../Kawabata_1981/PS/stress_stretch.txt");
data_3 = importdata("../Meunier_2008/PS/stress_stretch.txt");
P_exp_3 = data_3(:,1);
lambda_3 = data_3(:,2);
F_list_3 = zeros(3,3,length(lambda_3));
F_list_3(1,1,:) = lambda_3;
F_list_3(2,2,:) = 1.0;
F_list_3(3,3,:) = lambda_3.^(-1.0);

paras0 = [1.0, 50.0];
lb = [0.0, 0.0];
ub = [Inf, Inf];

objectiveFunction = @(paras) objective(paras, F_list_1, P_exp_1,...
                                              F_list_2, P_exp_2,...
                                              F_list_3, P_exp_3);
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
plot_results(paras,...
             F_list_1, P_exp_1,...
             F_list_2, P_exp_2,...
             F_list_3, P_exp_3);
print(gcf, '-djpeg', 'hyper_fit_AB.jpg');