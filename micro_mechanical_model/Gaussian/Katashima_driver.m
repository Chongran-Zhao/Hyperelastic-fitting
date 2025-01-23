clc; clear; close all;
addpath("src/")

% loading data of UT
lambda_ET = importdata('../Katashima_2012/ET/stretch.txt');
P1_ET = importdata('../Katashima_2012/ET/stress.txt');

F_ET_list = zeros(3,3,length(lambda_ET));
F_ET_list(1,1,:) = lambda_ET;
F_ET_list(2,2,:) = lambda_ET;
F_ET_list(3,3,:) = lambda_ET.^(-2.0);

paras0 = [0.05];
lb = [0.0];
ub = [Inf];
objectiveFunction = @(paras) objective(paras, F_ET_list, P1_ET);

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

% loading data of BT
lambda_BT_1 = importdata('../Katashima_2012/BT/stretch_1.txt');
relative_P1_BT_1 = importdata('../Katashima_2012/BT/stress_1.txt');
F_BT_list_1 = zeros(3,3,length(lambda_BT_1));
F_BT_list_1(1,1,:) = lambda_BT_1;
F_BT_list_1(2,2,:) = 0.5 + 0.5 .* lambda_BT_1;
F_BT_list_1(3,3,:) = 1.0 ./ ((0.5 + 0.5 .* lambda_BT_1) .* lambda_BT_1);

% loading data of BT
lambda_BT_2 = importdata('../Katashima_2012/BT/stretch_1_2.txt');
relative_P1_BT_2 = importdata('../Katashima_2012/BT/stress_2.txt');
F_BT_list_2 = zeros(3,3,length(lambda_BT_2));
F_BT_list_2(1,1,:) = lambda_BT_2;
F_BT_list_2(2,2,:) = 0.5 + 0.5 .* lambda_BT_2;
F_BT_list_2(3,3,:) = 1.0 ./ ((0.5 + 0.5 .* lambda_BT_2) .* lambda_BT_2);


% loading data of PS
lambda_PS_1 = importdata('../Katashima_2012/PS/stretch_1.txt');
relative_P1_PS_1 = importdata('../Katashima_2012/PS/stress_1.txt');
F_PS_list_1 = zeros(3,3,length(lambda_PS_1));
F_PS_list_1(1,1,:) = lambda_PS_1;
F_PS_list_1(2,2,:) = 1.0;
F_PS_list_1(3,3,:) = 1.0 ./ lambda_PS_1;

% loading data of PS
lambda_PS_2 = importdata('../Katashima_2012/PS/stretch_1_2.txt');
relative_P1_PS_2 = importdata('../Katashima_2012/PS/stress_2.txt');
F_PS_list_2 = zeros(3,3,length(lambda_PS_2));
F_PS_list_2(1,1,:) = lambda_PS_2;
F_PS_list_2(2,2,:) = 1.0;
F_PS_list_2(3,3,:) = 1.0 ./ lambda_PS_2;

plot_results_Katashima(paras,...
             F_ET_list, P1_ET, ...
             F_BT_list_1, relative_P1_BT_1, ...
             F_BT_list_2, relative_P1_BT_2, ...
             F_PS_list_1, relative_P1_PS_1, ...
             F_PS_list_2, relative_P1_PS_2);