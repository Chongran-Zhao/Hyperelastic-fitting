%% This code is to draw figures by parameters input
% Ogden Model with 6 parameters and GS Model with 6 parameters
% Input variables are Ogden_paras and GS_paras

clc;clear;
close all

% loading data
Treloar_UT_strain = importdata("./Treloar-UT/strain.txt");
Treloar_UT_stress = importdata("./Treloar-UT/stress.txt");

Treloar_ET_strain = importdata("./Treloar-ET/strain.txt");
Treloar_ET_stress = importdata("./Treloar-ET/stress.txt");

Treloar_PS_strain = importdata("./Treloar-PS/strain.txt");
Treloar_PS_stress = importdata("./Treloar-PS/stress.txt");

% Ogden_paras = [mu_1, alpha_1, mu_2, alpha_2, mu_3, alpha_3]
% Unit of mu: MPa
% Ogden_paras = [4.05965e-1, 1.819059, -1.879932e-3, -2.595254, 1.822613e-6, 8.189521];
Ogden_paras = [-0.00399450312646, -2.357101026245013, 0.000056818096938, 6.540469509568156, 0.516156037978041, 1.559863937925930];
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

% GS_paras = [n_1, m_1, mu_1, n_2, m_2, mu_2] 
% Unit of mu: MPa
GS_paras = [-0.961211379695161, -0.457753246616541, 0.353505442226954, 2.384084327573897, 3.989939972681667, 0.000013890096421];

% tool function for generalized strain
term1 = @(x, xdata) 2*x(3)*(xdata.^x(2) - xdata.^(-x(1))) .* ((x(2).*(xdata.^(x(2)-1)) + x(1).*(xdata.^(-x(1)-1)) )  / (x(2)+x(1)).^2);
term2 = @(x, xdata) 2*x(6)*(xdata.^x(5) - xdata.^(-x(4))) .* ((x(5).*(xdata.^(x(5)-1)) + x(4).*(xdata.^(-x(4)-1)) )  / (x(5)+x(4)).^2);

% P_11 of generalized strain
GS_UT = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-1.5)) .* ( term1(x, xdata.^(-0.5)) + term2(x, xdata.^(-0.5)) );
GS_ET = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-3.0)) .* ( term1(x, xdata.^(-2.0)) + term2(x, xdata.^(-2.0)) );
GS_PS = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-2.0)) .* ( term1(x, xdata.^(-1.0)) + term2(x, xdata.^(-1.0)) );

% Residual norm
GS_res = sum((GS_UT(GS_paras, Treloar_UT_strain)-Treloar_UT_stress) .* (GS_UT(GS_paras, Treloar_UT_strain)-Treloar_UT_stress) ./ Treloar_UT_stress)...
       + sum((GS_ET(GS_paras, Treloar_ET_strain)-Treloar_ET_stress) .* (GS_ET(GS_paras, Treloar_ET_strain)-Treloar_ET_stress) ./ Treloar_ET_stress)...
       + sum((GS_PS(GS_paras, Treloar_PS_strain)-Treloar_PS_stress) .* (GS_PS(GS_paras, Treloar_PS_strain)-Treloar_PS_stress) ./ Treloar_PS_stress);

disp(['GS residual norm = ' num2str(GS_res)]);

Ogden_res = sum((Ogden_UT(Ogden_paras, Treloar_UT_strain)-Treloar_UT_stress) .* (Ogden_UT(Ogden_paras, Treloar_UT_strain)-Treloar_UT_stress) ./ Treloar_UT_stress)...
          + sum((Ogden_ET(Ogden_paras, Treloar_ET_strain)-Treloar_ET_stress) .* (Ogden_ET(Ogden_paras, Treloar_ET_strain)-Treloar_ET_stress) ./ Treloar_ET_stress)...
          + sum((Ogden_PS(Ogden_paras, Treloar_PS_strain)-Treloar_PS_stress) .* (Ogden_PS(Ogden_paras, Treloar_PS_strain)-Treloar_PS_stress) ./ Treloar_PS_stress);

disp(['Ogden residual norm = ' num2str(Ogden_res)]);


% Prescribed stretch
x_UT = linspace(1.0, max(Treloar_UT_strain), 25);
x_ET = linspace(1.0, max(Treloar_ET_strain), 25);
x_PS = linspace(1.0, max(Treloar_PS_strain), 25);

% Plot GS
figure;
set(gcf, 'Position', [100, 100, 1000, 1000]); 
hold on;

plot(Treloar_UT_strain, Treloar_UT_stress, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(x_UT, GS_UT(GS_paras, x_UT), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');

plot(Treloar_ET_strain, Treloar_ET_stress, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(x_ET, GS_ET(GS_paras, x_ET), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-');

plot(Treloar_PS_strain, Treloar_PS_stress, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 12); 
plot(x_PS, GS_PS(GS_paras, x_PS), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '-');

hold off;

legend('UT Experimental', 'UT Generalized Strain Energy Model',...
       'ET Experimental', 'ET Generalized Strain Energy Model',...
       'PS Experimental', 'PS Generalized Strain Energy Model',...
       'Location', 'northwest', 'FontSize', 18);

title('Generalized Strain Energy Model Fitted Curves', 'FontSize', 20);
xlabel('Stretch', 'FontSize', 20);
ylabel('P_{11} Stress (MPa)', 'FontSize', 20);

saveas(gcf, 'GS.png');

% Plot Ogden
figure;
set(gcf, 'Position', [100, 100, 1000, 1000]); 
hold on;

plot(Treloar_UT_strain, Treloar_UT_stress, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(x_UT, Ogden_UT(Ogden_paras, x_UT), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');

plot(Treloar_ET_strain, Treloar_ET_stress, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(x_ET, Ogden_ET(Ogden_paras, x_ET), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-');

plot(Treloar_PS_strain, Treloar_PS_stress, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 12); 
plot(x_PS, Ogden_PS(Ogden_paras, x_PS), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '-');

hold off;

legend('UT Experimental', 'UT Ogden Model',...
       'ET Experimental', 'ET Ogden Model',...
       'PS Experimental', 'PS Ogden Model',...
       'Location', 'northwest', 'FontSize', 18);

title('Ogden Model Fitted Curves', 'FontSize', 20);
xlabel('Stretch', 'FontSize', 20);
ylabel('P_{11} Stress (MPa)', 'FontSize', 20);

saveas(gcf, 'Ogden.png');

