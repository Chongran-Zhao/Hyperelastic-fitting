clc;clear;
close all

% Prescribed stretch
UT_x = linspace(1, 8, 100);
ET_x = linspace(1, 5, 100);
PS_x = linspace(1, 5, 100);

% Ogden_paras = [mu_1, alpha_1, mu_2, alpha_2, mu_3, alpha_3]
% Unit of mu: Psi
Ogden_paras = [20, 1.8, -7, -2, 1.5, 7];
% Ogden_paras = [0.63, 1.3, -0.01, -2, 0.0012, 5];
mu = 0.5 * ( Ogden_paras(1)*Ogden_paras(2) + Ogden_paras(3)*Ogden_paras(4) + Ogden_paras(5)*Ogden_paras(6) );

% P_11 of Ogden Model
Ogden_UT = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-0.5 * x(2) - 1.0) ) ... 
+ x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-0.5 * x(4) - 1.0) ) ...
+ x(5) * ( xdata .^ (x(6) - 1.0) - xdata .^ (-0.5 * x(6) - 1.0) );

Ogden_ET = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-2.0 * x(2) - 1.0) ) ... 
+ x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-2.0 * x(4) - 1.0) ) ...
+ x(5) * ( xdata .^ (x(6) - 1.0) - xdata .^ (-2.0 * x(6) - 1.0) );

Ogden_PS = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-1.0 * x(2) - 1.0) ) ... 
+ x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-1.0 * x(4) - 1.0) ) ...
+ x(5) * ( xdata .^ (x(6) - 1.0) - xdata .^ (-1.0 * x(6) - 1.0) );

% Initial guess
% GS_paras = [n_1, m_1, mu_1, n_2, m_2, mu_2] 
% Unit of mu: Psi
GS_paras_0 = [1.0, 1.0, 25, 0.1, 1.0];

% tool function for generalized strain
term1 = @(x, xdata) 2*x(3)*(xdata.^x(2) - xdata.^(-x(1))) .* ((x(2).*(xdata.^(x(2)-1)) + x(1).*(xdata.^(-x(1)-1)) )  / (x(2)+x(1)).^2);
term2 = @(x, xdata) 2* (mu-x(3)) *(xdata.^x(5) - xdata.^(-x(4))) .* ((x(5).*(xdata.^(x(5)-1)) + x(4).*(xdata.^(-x(4)-1)) )  / (x(5)+x(4)).^2);

% P_11 of generalized strain
GS_UT = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-1.5)) .* ( term1(x, xdata.^(-0.5)) + term2(x, xdata.^(-0.5)) );
GS_ET = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-3.0)) .* ( term1(x, xdata.^(-2.0)) + term2(x, xdata.^(-2.0)) );
GS_PS = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-2.0)) .* ( term1(x, xdata.^(-1.0)) + term2(x, xdata.^(-1.0)) );

objectiveFunction = @(x) sum((GS_UT(x, UT_x) - Ogden_UT(Ogden_paras, UT_x)).^2) ./ length(UT_x) +... 
                         sum((GS_ET(x, ET_x) - Ogden_ET(Ogden_paras, ET_x)).^2) ./ length(ET_x) +...
                         sum((GS_PS(x, PS_x) - Ogden_PS(Ogden_paras, PS_x)).^2) ./ length(PS_x);

lb = [-Inf, -Inf, 0, -Inf, -Inf];
ub = [Inf, Inf, mu, Inf, Inf];

options = optimoptions('lsqnonlin', 'Algorithm', 'interior-point', 'MaxIterations', 5000);
[GS_paras, ~] = lsqnonlin( objectiveFunction, GS_paras_0, lb, ub, options );

disp(['n_1 = ' num2str(GS_paras(1)) ' m_1 = ' num2str(GS_paras(2)) ' mu_1 = ' num2str(GS_paras(3))]);
disp(['n_2 = ' num2str(GS_paras(4)) ' m_2 = ' num2str(GS_paras(5)) ' mu_2 = ' num2str(mu-GS_paras(3))]);
resnorm = objectiveFunction(GS_paras);
disp(['residual norm = ' num2str(resnorm)] );

% figures
figure;
hold on;

plot(UT_x, GS_UT(GS_paras, UT_x), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2); % 深灰色
plot(UT_x, Ogden_UT(Ogden_paras, UT_x), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2, 'LineStyle', '--'); % 浅灰色

plot(ET_x, GS_ET(GS_paras, ET_x), 'Color', [0.85, 0.33, 0], 'LineWidth', 2); % 橙色
plot(ET_x, Ogden_ET(Ogden_paras, ET_x), 'Color', [0.93, 0.69, 0.13], 'LineWidth', 2, 'LineStyle', '--'); % 金色

plot(PS_x, GS_PS(GS_paras, PS_x), 'Color', [0, 0.5, 0.5], 'LineWidth', 2); % 青色
plot(PS_x, Ogden_PS(Ogden_paras, PS_x), 'Color', [0.27, 0.51, 0.71], 'LineWidth', 2, 'LineStyle', '--'); % 天蓝色

hold off;
grid off;

legend('UT fitted by Generalized strain', 'UT fitted by Ogden Model',...
    'ET fitted by Generalized strain', 'ET fitted by Ogden Model',...
    'PS fitted by Generalized strain', 'PS fitted by Ogden Model',...
    'Location', 'NorthWest', 'FontSize', 14);

title('Data and Fitted Curves', 'FontSize', 14);
xlabel('Stretch', 'FontSize', 12);
ylabel('P_{11} of 1st PK stress', 'FontSize', 12);
