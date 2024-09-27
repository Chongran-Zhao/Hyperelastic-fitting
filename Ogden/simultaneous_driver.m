clc; clear; close all;
addpath("src/")

% loading data of UT
data_1 = importdata("../Treloar_1944/UT/stress_stretch.txt");
% data_1 = importdata("../Kawabata_1981/UT/stress_stretch.txt");
% data_1 = importdata("../Meunier_2008/UT/stress_stretch.txt");
P_exp_1 = data_1(:,1);
lambda_1 = data_1(:,2);
Ft_1 = zeros(3,3,length(lambda_1));
Ft_1(1,1,:) = lambda_1;
Ft_1(2,2,:) = lambda_1.^(-0.5);
Ft_1(3,3,:) = lambda_1.^(-0.5);

% loading data of ET
data_2 = importdata("../Treloar_1944/ET/stress_stretch.txt");
% data_2 = importdata("../Kawabata_1981/ET/stress_stretch.txt");
% data_2 = importdata("../Meunier_2008/ET/stress_stretch.txt");
P_exp_2 = data_2(:,1);
lambda_2 = data_2(:,2);
Ft_2 = zeros(3,3,length(lambda_2));
Ft_2(1,1,:) = lambda_2;
Ft_2(2,2,:) = lambda_2;
Ft_2(3,3,:) = lambda_2.^(-2.0);

% loading data of PS
data_3 = importdata("../Treloar_1944/PS/stress_stretch.txt");
% data_3 = importdata("../Kawabata_1981/PS/stress_stretch.txt");
% data_3 = importdata("../Meunier_2008/PS/stress_stretch.txt");
P_exp_3 = data_3(:,1);
lambda_3 = data_3(:,2);
Ft_3 = zeros(3,3,length(lambda_3));
Ft_3(1,1,:) = lambda_3;
Ft_3(2,2,:) = 1.0;
Ft_3(3,3,:) = lambda_3.^(-1.0);

num_terms = 3;
paras0 = [1.0, 1.0, -1.0, ...
         1.0, 1.0, -1.0 ];
lb = [0.0, 0.0, -Inf, ...
      0.0, 0.0, -Inf ];
ub = [Inf, Inf, 0.0, ...
      Inf, Inf, 0.0 ];

objectiveFunction = @(paras) objective(paras, num_terms, Ft_1, P_exp_1,...
                                                         Ft_2, P_exp_2,...
                                                         Ft_3, P_exp_3);
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 40000,  ...
    'TolFun', 1.0e-10, ...
    'TolX', 1.0e-10, ...
    'OptimalityTolerance', 1.0e-10, ...
    'StepTolerance', 1.0e-10, ...
    'Display', 'iter');

[paras, ~] = lsqnonlin( objectiveFunction, paras0, [], [], options);
plot_results(paras, num_terms,...
             Ft_1, P_exp_1,...
             Ft_2, P_exp_2,...
             Ft_3, P_exp_3);
print(gcf, '-djpeg', 'hyper_fit_Ogden.jpg');