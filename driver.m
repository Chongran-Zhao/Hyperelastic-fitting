clc;clear;
addpath('src/');
% loading data
lambda1_UT = importdata("./Treloar-UT/strain.txt");
P1_exp_UT = importdata("./Treloar-UT/stress.txt");

lambda1_ET = importdata("./Treloar-ET/strain.txt");
P1_exp_ET = importdata("./Treloar-ET/stress.txt");

lambda1_PS = importdata("./Treloar-PS/strain.txt");
P1_exp_PS = importdata("./Treloar-PS/stress.txt");

input_num_paras = 6;
paras0 = ones(1,input_num_paras);

Model = Ogden(input_num_paras, lambda1_UT, lambda1_ET, lambda1_PS, P1_exp_UT, P1_exp_ET, P1_exp_PS);

objectiveFunction = @(paras) Model.objective(paras);

% options = optimoptions('lsqnonlin', 'Algorithm', 'interior-point', 'MaxIterations', 5000);
% [paras_opt, res] = lsqnonlin(objectiveFunction, paras0, [], [], options);