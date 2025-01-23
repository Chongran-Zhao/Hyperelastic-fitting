clc; clear; close all;
addpath("src/")

% loading data of UT
lambda_UT = importdata('../Jones_1975/UT/stretch.txt');
P1_UT = importdata('../Jones_1975/UT/stress.txt');

F_UT_list = zeros(3,3,length(lambda_UT));
F_UT_list(1,1,:) = lambda_UT;
F_UT_list(2,2,:) = lambda_UT.^(-0.5);
F_UT_list(3,3,:) = lambda_UT.^(-0.5);

paras0 = [1.0, 25.0];
lb = [0.0, 0.0];
ub = [Inf, Inf];

objectiveFunction = @(paras) objective(paras, F_UT_list, P1_UT);

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
lambda_BT_1 = importdata('../Jones_1975/Biaxial_tension/lambda_1/stretch.txt');
relative_P1_BT_1 = importdata('../Jones_1975/Biaxial_tension/lambda_1/relative_stress.txt');
F_BT_list_1 = zeros(3,3,length(lambda_BT_1));
F_BT_list_1(1,1,:) = lambda_BT_1;
F_BT_list_1(2,2,:) = ones(length(lambda_BT_1), 1);
F_BT_list_1(3,3,:) = 1.0 ./ (1.0 .* lambda_BT_1);

% loading data of BT
lambda_BT_2 = importdata('../Jones_1975/Biaxial_tension/lambda_1d502/stretch.txt');
relative_P1_BT_2 = importdata('../Jones_1975/Biaxial_tension/lambda_1d502/relative_stress.txt');
F_BT_list_2 = zeros(3,3,length(lambda_BT_2));
F_BT_list_2(1,1,:) = lambda_BT_2;
F_BT_list_2(2,2,:) = 1.502 .* ones(length(lambda_BT_2), 1);
F_BT_list_2(3,3,:) = 1.0 ./ (1.502 .* lambda_BT_2);

% loading data of BT
lambda_BT_3 = importdata('../Jones_1975/Biaxial_tension/lambda_1d984/stretch.txt');
relative_P1_BT_3 = importdata('../Jones_1975/Biaxial_tension/lambda_1d984/relative_stress.txt');
F_BT_list_3 = zeros(3,3,length(lambda_BT_3));
F_BT_list_3(1,1,:) = lambda_BT_3;
F_BT_list_3(2,2,:) = 1.984 .* ones(length(lambda_BT_3), 1);
F_BT_list_3(3,3,:) = 1.0 ./ (1.984 .* lambda_BT_3);

% loading data of BT
lambda_BT_4 = importdata('../Jones_1975/Biaxial_tension/lambda_2d295/stretch.txt');
relative_P1_BT_4 = importdata('../Jones_1975/Biaxial_tension/lambda_2d295/relative_stress.txt');
F_BT_list_4 = zeros(3,3,length(lambda_BT_4));
F_BT_list_4(1,1,:) = lambda_BT_4;
F_BT_list_4(2,2,:) = 2.295 .* ones(length(lambda_BT_4), 1);
F_BT_list_4(3,3,:) = 1.0 ./ (2.295 .* lambda_BT_4);

% loading data of BT
lambda_BT_5 = importdata('../Jones_1975/Biaxial_tension/lambda_2d623/stretch.txt');
relative_P1_BT_5 = importdata('../Jones_1975/Biaxial_tension/lambda_2d623/relative_stress.txt');
F_BT_list_5 = zeros(3,3,length(lambda_BT_5));
F_BT_list_5(1,1,:) = lambda_BT_5;
F_BT_list_5(2,2,:) = 2.623 .* ones(length(lambda_BT_5), 1);
F_BT_list_5(3,3,:) = 1.0 ./ (2.623 .* lambda_BT_5);

plot_results_Jones(paras,...
             F_UT_list, P1_UT, ...
             F_BT_list_1, relative_P1_BT_1, ...
             F_BT_list_2, relative_P1_BT_2, ...
             F_BT_list_3, relative_P1_BT_3, ...
             F_BT_list_4, relative_P1_BT_4, ...
             F_BT_list_5, relative_P1_BT_5);
% print(gcf, '-djpeg', 'hyper_fit_Zhan.jpg');