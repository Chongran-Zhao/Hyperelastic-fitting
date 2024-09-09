clc; clear; close all;
addpath("src/")

% loading data of UT
lambda_1 = importdata("../exp-data/Treloar-UT/strain.txt");
P_exp_1 = importdata("../exp-data/Treloar-UT/stress.txt");
F_list_1 = zeros(3,3,length(lambda_1));
F_list_1(1,1,:) = lambda_1;
F_list_1(2,2,:) = lambda_1.^(-0.5);
F_list_1(3,3,:) = lambda_1.^(-0.5);

% loading data of ET
lambda_2 = importdata("../exp-data/Treloar-ET/strain.txt");
P_exp_2 = importdata("../exp-data/Treloar-ET/stress.txt");
F_list_2 = zeros(3,3,length(lambda_2));
F_list_2(1,1,:) = lambda_2;
F_list_2(2,2,:) = lambda_2;
F_list_2(3,3,:) = lambda_2.^(-2.0);

% loading data of UT
lambda_3 = importdata("../exp-data/Treloar-PS/strain.txt");
P_exp_3 = importdata("../exp-data/Treloar-PS/stress.txt");
F_list_3 = zeros(3,3,length(lambda_3));
F_list_3(1,1,:) = lambda_3;
F_list_3(2,2,:) = 1.0;
F_list_3(3,3,:) = lambda_3.^(-1.0);

xi_0 = [1.0, 1.0, 1.0];
objectiveFunction = @(xi) objective(xi, F_list_1, P_exp_1,...
                                        F_list_2, P_exp_2,...
                                        F_list_3, P_exp_3);
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxIterations', 1000, ...
    'Display', 'iter');

[xi, ~] = lsqnonlin( objectiveFunction, xi_0, [], [], options);
plot_results(xi,...
             F_list_1, P_exp_1,...
             F_list_2, P_exp_2,...
             F_list_3, P_exp_3);
print(gcf, '-djpeg', 'hyper_fit_Yeoh.jpg');