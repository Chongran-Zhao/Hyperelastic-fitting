clc; clear; close all;
addpath("src/")

% loading data of UT
lambda_1 = importdata("../experimental-data/Treloar-UT/strain.txt");
P_exp_1 = importdata("../experimental-data/Treloar-UT/stress.txt");
Ft_1 = zeros(3,3,length(lambda_1));
Ft_1(1,1,:) = lambda_1;
Ft_1(2,2,:) = lambda_1.^(-0.5);
Ft_1(3,3,:) = lambda_1.^(-0.5);

% loading data of ET
lambda_2 = importdata("../experimental-data/Treloar-ET/strain.txt");
P_exp_2 = importdata("../experimental-data/Treloar-ET/stress.txt");
Ft_2 = zeros(3,3,length(lambda_2));
Ft_2(1,1,:) = lambda_2;
Ft_2(2,2,:) = lambda_2;
Ft_2(3,3,:) = lambda_2.^(-2.0);

% loading data of UT
lambda_3 = importdata("../experimental-data/Treloar-PS/strain.txt");
P_exp_3 = importdata("../experimental-data/Treloar-PS/stress.txt");
Ft_3 = zeros(3,3,length(lambda_3));
Ft_3(1,1,:) = lambda_3;
Ft_3(2,2,:) = 1.0;
Ft_3(3,3,:) = lambda_3.^(-1.0);

xi0 = [1.0, 1.0, 1.0];
objectiveFunction = @(xi) multi_objective(xi, Ft_1, P_exp_1,...
                                      Ft_2, P_exp_2,...
                                      Ft_3, P_exp_3);
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxIterations', 1000, ...
    'Display', 'iter');

[xi, ~] = lsqnonlin( objectiveFunction, xi0, [], [], options);
plot_results(xi,...
    Ft_1, P_exp_1,...
    Ft_2, P_exp_2,...
    Ft_3, P_exp_3);
print(gcf, '-djpeg', 'fig_sim.jpg');