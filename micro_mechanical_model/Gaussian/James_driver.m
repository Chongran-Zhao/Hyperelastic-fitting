clc; clear; close all;
addpath("src/")

% loading data of UT
lambda_UT = importdata('../James_1975/UT/stretch.txt');
P1_UT = importdata('../James_1975/UT/stress.txt');

F_UT_list = zeros(3,3,length(lambda_UT));
F_UT_list(1,1,:) = lambda_UT;
F_UT_list(2,2,:) = lambda_UT.^(-0.5);
F_UT_list(3,3,:) = lambda_UT.^(-0.5);

paras0 = [0.05];
lb = [0.0];
ub = [Inf];
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
lambda_BT_1 = importdata('../James_1975/BT/lambda_3d5/stretch_1.txt');
relative_P1_BT_1 = importdata('../James_1975/BT/lambda_3d5/stress_1.txt');
F_BT_list_1 = zeros(3,3,length(lambda_BT_1));
F_BT_list_1(1,1,:) = lambda_BT_1;
F_BT_list_1(2,2,:) = 3.5 .* ones(length(lambda_BT_1), 1);
F_BT_list_1(3,3,:) = 1.0 ./ (3.5 .* lambda_BT_1);

% loading data of BT
lambda_BT_2 = importdata('../James_1975/BT/lambda_3d0/stretch_1.txt');
relative_P1_BT_2 = importdata('../James_1975/BT/lambda_3d0/stress_1.txt');
F_BT_list_2 = zeros(3,3,length(lambda_BT_2));
F_BT_list_2(1,1,:) = lambda_BT_2;
F_BT_list_2(2,2,:) = 3.0 .* ones(length(lambda_BT_2), 1);
F_BT_list_2(3,3,:) = 1.0 ./ (3.0 .* lambda_BT_2);

% loading data of BT
lambda_BT_3 = importdata('../James_1975/BT/lambda_2d5/stretch_1.txt');
relative_P1_BT_3 = importdata('../James_1975/BT/lambda_2d5/stress_1.txt');
F_BT_list_3 = zeros(3,3,length(lambda_BT_3));
F_BT_list_3(1,1,:) = lambda_BT_3;
F_BT_list_3(2,2,:) = 2.5 .* ones(length(lambda_BT_3), 1);
F_BT_list_3(3,3,:) = 1.0 ./ (2.5 .* lambda_BT_3);

% loading data of BT
lambda_BT_4 = importdata('../James_1975/BT/lambda_2d0/stretch_1.txt');
relative_P1_BT_4 = importdata('../James_1975/BT/lambda_2d0/stress_1.txt');
F_BT_list_4 = zeros(3,3,length(lambda_BT_4));
F_BT_list_4(1,1,:) = lambda_BT_4;
F_BT_list_4(2,2,:) = 2.0 .* ones(length(lambda_BT_4), 1);
F_BT_list_4(3,3,:) = 1.0 ./ (2.0 .* lambda_BT_4);

lambda_BT_5 = importdata('../James_1975/BT/lambda_1d7/stretch_1.txt');
relative_P1_BT_5 = importdata('../James_1975/BT/lambda_1d7/stress_1.txt');
F_BT_list_5 = zeros(3,3,length(lambda_BT_5));
F_BT_list_5(1,1,:) = lambda_BT_5;
F_BT_list_5(2,2,:) = 1.7 .* ones(length(lambda_BT_5), 1);
F_BT_list_5(3,3,:) = 1.0 ./ (1.7 .* lambda_BT_5);

lambda_BT_6 = importdata('../James_1975/BT/lambda_1d5/stretch_1.txt');
relative_P1_BT_6 = importdata('../James_1975/BT/lambda_1d5/stress_1.txt');
F_BT_list_6 = zeros(3,3,length(lambda_BT_6));
F_BT_list_6(1,1,:) = lambda_BT_6;
F_BT_list_6(2,2,:) = 1.5 .* ones(length(lambda_BT_6), 1);
F_BT_list_6(3,3,:) = 1.0 ./ (1.5 .* lambda_BT_6);

lambda_BT_7 = importdata('../James_1975/BT/lambda_1d3/stretch_1.txt');
relative_P1_BT_7 = importdata('../James_1975/BT/lambda_1d3/stress_1.txt');
F_BT_list_7 = zeros(3,3,length(lambda_BT_7));
F_BT_list_7(1,1,:) = lambda_BT_7;
F_BT_list_7(2,2,:) = 1.3 .* ones(length(lambda_BT_7), 1);
F_BT_list_7(3,3,:) = 1.0 ./ (1.3 .* lambda_BT_7);

lambda_BT_8 = importdata('../James_1975/BT/lambda_3d5/stretch_1_2.txt');
relative_P1_BT_8 = importdata('../James_1975/BT/lambda_3d5/stress_2.txt');
F_BT_list_8 = zeros(3,3,length(lambda_BT_8));
F_BT_list_8(1,1,:) = lambda_BT_8;
F_BT_list_8(2,2,:) = 3.5 .* ones(length(lambda_BT_8), 1);
F_BT_list_8(3,3,:) = 1.0 ./ (3.5 .* lambda_BT_8);

% loading data of BT
lambda_BT_9 = importdata('../James_1975/BT/lambda_3d0/stretch_1_2.txt');
relative_P1_BT_9 = importdata('../James_1975/BT/lambda_3d0/stress_2.txt');
F_BT_list_9 = zeros(3,3,length(lambda_BT_9));
F_BT_list_9(1,1,:) = lambda_BT_9;
F_BT_list_9(2,2,:) = 3.0 .* ones(length(lambda_BT_9), 1);
F_BT_list_9(3,3,:) = 1.0 ./ (3.0 .* lambda_BT_9);

% loading data of BT
lambda_BT_10 = importdata('../James_1975/BT/lambda_2d5/stretch_1_2.txt');
relative_P1_BT_10 = importdata('../James_1975/BT/lambda_2d5/stress_2.txt');
F_BT_list_10 = zeros(3,3,length(lambda_BT_10));
F_BT_list_10(1,1,:) = lambda_BT_10;
F_BT_list_10(2,2,:) = 2.5 .* ones(length(lambda_BT_10), 1);
F_BT_list_10(3,3,:) = 1.0 ./ (2.5 .* lambda_BT_10);

% loading data of BT
lambda_BT_11 = importdata('../James_1975/BT/lambda_2d0/stretch_1_2.txt');
relative_P1_BT_11 = importdata('../James_1975/BT/lambda_2d0/stress_2.txt');
F_BT_list_11 = zeros(3,3,length(lambda_BT_11));
F_BT_list_11(1,1,:) = lambda_BT_11;
F_BT_list_11(2,2,:) = 2.0 .* ones(length(lambda_BT_11), 1);
F_BT_list_11(3,3,:) = 1.0 ./ (2.0 .* lambda_BT_11);

% loading data of BT
lambda_BT_12 = importdata('../James_1975/BT/lambda_1d7/stretch_1_2.txt');
relative_P1_BT_12 = importdata('../James_1975/BT/lambda_1d7/stress_2.txt');
F_BT_list_12 = zeros(3,3,length(lambda_BT_12));
F_BT_list_12(1,1,:) = lambda_BT_12;
F_BT_list_12(2,2,:) = 1.7 .* ones(length(lambda_BT_12), 1);
F_BT_list_12(3,3,:) = 1.0 ./ (1.7 .* lambda_BT_12);

% loading data of BT
lambda_BT_13 = importdata('../James_1975/BT/lambda_1d5/stretch_1_2.txt');
relative_P1_BT_13 = importdata('../James_1975/BT/lambda_1d5/stress_2.txt');
F_BT_list_13 = zeros(3,3,length(lambda_BT_13));
F_BT_list_13(1,1,:) = lambda_BT_13;
F_BT_list_13(2,2,:) = 1.5 .* ones(length(lambda_BT_13), 1);
F_BT_list_13(3,3,:) = 1.0 ./ (1.5 .* lambda_BT_13);

% loading data of BT
lambda_BT_14 = importdata('../James_1975/BT/lambda_1d3/stretch_1_2.txt');
relative_P1_BT_14 = importdata('../James_1975/BT/lambda_1d3/stress_2.txt');
F_BT_list_14 = zeros(3,3,length(lambda_BT_14));
F_BT_list_14(1,1,:) = lambda_BT_14;
F_BT_list_14(2,2,:) = 1.3 .* ones(length(lambda_BT_14), 1);
F_BT_list_14(3,3,:) = 1.0 ./ (1.3 .* lambda_BT_14);

plot_results_James(paras,...
             F_UT_list, P1_UT, ...
             F_BT_list_1, relative_P1_BT_1, ...
             F_BT_list_2, relative_P1_BT_2, ...
             F_BT_list_3, relative_P1_BT_3, ...
             F_BT_list_4, relative_P1_BT_4, ...
             F_BT_list_5, relative_P1_BT_5, ...
             F_BT_list_6, relative_P1_BT_6, ...
             F_BT_list_7, relative_P1_BT_7, ...
             F_BT_list_8, relative_P1_BT_8, ...
             F_BT_list_9, relative_P1_BT_9, ...
             F_BT_list_10, relative_P1_BT_10, ...
             F_BT_list_11, relative_P1_BT_11, ...
             F_BT_list_12, relative_P1_BT_12, ...
             F_BT_list_13, relative_P1_BT_13, ...
             F_BT_list_14, relative_P1_BT_14);