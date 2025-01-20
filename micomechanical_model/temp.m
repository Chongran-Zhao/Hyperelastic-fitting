clc; clear; close all;

paras = [20.0, 100.0];
C = eye(3);  % 单位矩阵

result = non_Gaussian_energy(paras, C);
disp(['数值积分结果: ', num2str(result)]);

% 理论解
mu = paras(1);
N = paras(2);
beta = inv_Langevin(1 / sqrt(N));  % 假设 beta 是常数
theoretical_result = 4 * pi * mu * N * (beta / sqrt(N) + log(beta / sinh(beta)));
disp(['理论解: ', num2str(theoretical_result)]);