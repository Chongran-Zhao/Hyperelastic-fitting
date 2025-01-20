clc; clear; close all;
addpath("src/")

% loading data of UT
lambda_UT = importdata('../Kawamura_2001/UT/stretch.txt');
P1_UT = importdata('../Kawamura_2001/UT/stress.txt');

F_UT_list = zeros(3,3,length(lambda_UT));
F_UT_list(1,1,:) = lambda_UT;
F_UT_list(2,2,:) = lambda_UT.^(-0.5);
F_UT_list(3,3,:) = lambda_UT.^(-0.5);

% loading data of UE
lambda_UE = importdata('../Kawamura_2001/UE/stretch.txt');
P1_UE = importdata('../Kawamura_2001/UE/stress.txt');

F_UE_list = zeros(3,3,length(lambda_UE));
F_UE_list(1,1,:) = lambda_UE;
F_UE_list(2,2,:) = lambda_UE.^(-0.5);
F_UE_list(3,3,:) = lambda_UE.^(-0.5);

% loading data of ET
lambda_ET = importdata('../Kawamura_2001/ET/stretch.txt');
P1_ET = importdata('../Kawamura_2001/ET/stress.txt');

F_ET_list = zeros(3,3,length(lambda_ET));
F_ET_list(1,1,:) = lambda_ET;
F_ET_list(2,2,:) = lambda_ET;
F_ET_list(3,3,:) = lambda_ET.^(-2.0);

% paras0 = [0.1];
% lb = [0.0];
% ub = [Inf];
paras0 = [1.1, 100];
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
lambda_BT_1 = importdata('../Kawamura_2001/BT/lambda_1d7/stretch_1.txt');
relative_P1_BT_1 = importdata('../Kawamura_2001/BT/lambda_1d7/stress_1.txt');
F_BT_list_1 = zeros(3,3,length(lambda_BT_1));
F_BT_list_1(1,1,:) = lambda_BT_1;
F_BT_list_1(2,2,:) = 1.7 .* ones(length(lambda_BT_1), 1);
F_BT_list_1(3,3,:) = 1.0 ./ (1.7 .* lambda_BT_1);

% loading data of BT
lambda_BT_2 = importdata('../Kawamura_2001/BT/lambda_1d5/stretch_1.txt');
relative_P1_BT_2 = importdata('../Kawamura_2001/BT/lambda_1d5/stress_1.txt');
F_BT_list_2 = zeros(3,3,length(lambda_BT_2));
F_BT_list_2(1,1,:) = lambda_BT_2;
F_BT_list_2(2,2,:) = 1.5 .* ones(length(lambda_BT_2), 1);
F_BT_list_2(3,3,:) = 1.0 ./ (1.5 .* lambda_BT_2);

% loading data of BT
lambda_BT_3 = importdata('../Kawamura_2001/BT/lambda_1d3/stretch_1.txt');
relative_P1_BT_3 = importdata('../Kawamura_2001/BT/lambda_1d3/stress_1.txt');
F_BT_list_3 = zeros(3,3,length(lambda_BT_3));
F_BT_list_3(1,1,:) = lambda_BT_3;
F_BT_list_3(2,2,:) = 1.3 .* ones(length(lambda_BT_3), 1);
F_BT_list_3(3,3,:) = 1.0 ./ (1.3 .* lambda_BT_3);

% loading data of BT
lambda_BT_4 = importdata('../Kawamura_2001/BT/lambda_1d1/stretch_1.txt');
relative_P1_BT_4 = importdata('../Kawamura_2001/BT/lambda_1d1/stress_1.txt');
F_BT_list_4 = zeros(3,3,length(lambda_BT_4));
F_BT_list_4(1,1,:) = lambda_BT_4;
F_BT_list_4(2,2,:) = 1.1 .* ones(length(lambda_BT_4), 1);
F_BT_list_4(3,3,:) = 1.0 ./ (1.1 .* lambda_BT_4);

lambda_BT_5 = importdata('../Kawamura_2001/BT/lambda_1d7/stretch_1_2.txt');
relative_P1_BT_5 = importdata('../Kawamura_2001/BT/lambda_1d7/stress_2.txt');
F_BT_list_5 = zeros(3,3,length(lambda_BT_5));
F_BT_list_5(1,1,:) = lambda_BT_5;
F_BT_list_5(2,2,:) = 1.7 .* ones(length(lambda_BT_5), 1);
F_BT_list_5(3,3,:) = 1.0 ./ (1.7 .* lambda_BT_5);

lambda_BT_6 = importdata('../Kawamura_2001/BT/lambda_1d5/stretch_1_2.txt');
relative_P1_BT_6 = importdata('../Kawamura_2001/BT/lambda_1d5/stress_2.txt');
F_BT_list_6 = zeros(3,3,length(lambda_BT_6));
F_BT_list_6(1,1,:) = lambda_BT_6;
F_BT_list_6(2,2,:) = 1.5 .* ones(length(lambda_BT_6), 1);
F_BT_list_6(3,3,:) = 1.0 ./ (1.5 .* lambda_BT_6);

lambda_BT_7 = importdata('../Kawamura_2001/BT/lambda_1d3/stretch_1_2.txt');
relative_P1_BT_7 = importdata('../Kawamura_2001/BT/lambda_1d3/stress_2.txt');
F_BT_list_7 = zeros(3,3,length(lambda_BT_7));
F_BT_list_7(1,1,:) = lambda_BT_7;
F_BT_list_7(2,2,:) = 1.3 .* ones(length(lambda_BT_7), 1);
F_BT_list_7(3,3,:) = 1.0 ./ (1.3 .* lambda_BT_7);

lambda_BT_8 = importdata('../Kawamura_2001/BT/lambda_1d1/stretch_1_2.txt');
relative_P1_BT_8 = importdata('../Kawamura_2001/BT/lambda_1d1/stress_2.txt');
F_BT_list_8 = zeros(3,3,length(lambda_BT_8));
F_BT_list_8(1,1,:) = lambda_BT_8;
F_BT_list_8(2,2,:) = 1.1 .* ones(length(lambda_BT_8), 1);
F_BT_list_8(3,3,:) = 1.0 ./ (1.1 .* lambda_BT_8);

plot_results_Kawamura(paras,...
             F_UT_list, P1_UT, ...
             F_UE_list, P1_UE, ...
             F_ET_list, P1_ET, ...
             F_BT_list_1, relative_P1_BT_1, ...
             F_BT_list_2, relative_P1_BT_2, ...
             F_BT_list_3, relative_P1_BT_3, ...
             F_BT_list_4, relative_P1_BT_4, ...
             F_BT_list_5, relative_P1_BT_5, ...
             F_BT_list_6, relative_P1_BT_6, ...
             F_BT_list_7, relative_P1_BT_7, ...
             F_BT_list_8, relative_P1_BT_8);