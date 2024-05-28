clc; clear; close all;

% loading data
Treloar_UT_strain = importdata("./Treloar-UT/strain.txt");
Treloar_UT_stress = importdata("./Treloar-UT/stress.txt");

Treloar_ET_strain = importdata("./Treloar-ET/strain.txt");
Treloar_ET_stress = importdata("./Treloar-ET/stress.txt");

Treloar_PS_strain = importdata("./Treloar-PS/strain.txt");
Treloar_PS_stress = importdata("./Treloar-PS/stress.txt");
Model_name = 'MR Model';

[paras, UT, ET, PS, resnorm] = curve_fitting(Model_name, ...
                                    Treloar_UT_strain, Treloar_UT_stress, ...
                                    Treloar_ET_strain, Treloar_ET_stress, ...
                                    Treloar_PS_strain, Treloar_PS_stress);

UT_x = linspace(1.0, max(Treloar_UT_strain), 25);
ET_x = linspace(1.0, max(Treloar_ET_strain), 25);
PS_x = linspace(1.0, max(Treloar_PS_strain), 25);

figure;
hold on;

plot(Treloar_UT_strain, Treloar_UT_stress, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10); 
plot(UT_x, UT(paras, UT_x), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-');

plot(Treloar_ET_strain, Treloar_ET_stress, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10); 
plot(ET_x, ET(paras, ET_x), 'Color', 'r', 'LineWidth', 2);

plot(Treloar_PS_strain, Treloar_PS_stress, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10); 
plot(PS_x, PS(paras, PS_x), 'Color', 'b', 'LineWidth', 2);

format long
disp('Fitted Parameters:');

for i = 1:numel(paras)
    fprintf('para %d: %.6e\n', i, paras(i));
end
ylim( [0, 6] );
xlim( [1, 8] );


hXLabel = xlabel('Stretch', 'interpreter', 'latex');
hYLabel = ylabel('Nominal stress (MPa)', 'interpreter', 'latex');

set( gca, 'Box', 'on', 'TickDir'     , 'out', ...
    'TickLength'  , [.02 .02], ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'  , ...
    'YGrid'       , 'on' , ...
    'XGrid'       , 'on' , ...
    'XColor'      , [0 0 0 ], ...
    'YColor'      , [0 0 0 ], ...
    'LineWidth'   , 2 );
set(gca,'FontSize', 25,'fontWeight','bold');
set([hXLabel, hYLabel], 'FontName', 'Helvetica', 'FontSize', 30, 'FontWeight', 'bold');

l = legend( 'UT Exp', 'ET Exp', 'PS Exp', ...
            'UT Sim', 'ET Sim', 'PS Sim' );
set( l, 'interpreter','latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', 'Orientation', ...
     'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 3 );
l.ItemTokenSize=[45,45];
text(5.2, 0.5, ['$\chi^2$ = ', num2str(resnorm)], 'Interpreter', 'latex', 'FontSize', 40);
X = 42.0;
Y = X * 1.2;
xMargin = 3;
yMargin = 3;
xSize = X - 2 * xMargin;
ySize = Y - 2 * yMargin;
set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[X Y]);
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize]);
set(gcf, 'PaperOrientation','portrait');

saveas(gcf, [Model_name, '.png']);
% Curve fitting function
function [paras, UT, ET, PS, resnorm] = curve_fitting(Model_name, ...
                                             UT_strain, UT_stress, ...
                                             ET_strain, ET_stress, ...
                                             PS_strain, PS_stress)
switch Model_name
    case 'Ogden Model'
        [paras_0, lb, ub, UT, ET, PS] = Ogden_Model_Init();
    case 'Ogden4 Model'
        [paras_0, lb, ub, UT, ET, PS] = Ogden4_Model_Init();
    case 'CR Model'
        [paras_0, lb, ub, UT, ET, PS] = CR_Model_Init();
    case 'CR4 Model'
        [paras_0, lb, ub, UT, ET, PS] = CR4_Model_Init();
    case 'AB Model'
        [paras_0, lb, ub, UT, ET, PS] = AB_Model_Init();
    case 'MR Model'
        [paras_0, lb, ub, UT, ET, PS] = MR_Model_Init();
    otherwise
        error('ERROR: WRONG INPUT MODEL_NAME!');
end

objectiveFunction = @(xi) objective(xi, UT_strain, UT_stress, ET_strain, ET_stress, PS_strain, PS_stress, UT, ET, PS);
options = optimoptions('lsqnonlin', 'Algorithm', 'interior-point', 'MaxIterations', 5000);
[paras, ~] = lsqnonlin( objectiveFunction, paras_0, lb, ub, options);
resnorm = res(paras, UT_strain, UT_stress, ET_strain, ET_stress, PS_strain, PS_stress, UT, ET, PS);
disp(['Residual norm = ' num2str(resnorm)]);
value = objective(paras, UT_strain, UT_stress, ET_strain, ET_stress, PS_strain, PS_stress, UT, ET, PS);
disp(['Objective value = ' num2str(value)]);
end

% Objective function
function res = objective(xi, UT_strain, UT_stress, ET_strain, ET_stress, PS_strain, PS_stress, UT, ET, PS)

% mean square difference
% res = sum((UT(xi, UT_strain) - UT_stress).^2) ./ length(UT_strain)+ ...
%       sum((ET(xi, ET_strain) - ET_stress).^2) ./ length(ET_strain)+ ...
%       sum((PS(xi, PS_strain) - PS_stress).^2) ./ length(PS_strain);

% Relative Error (RE)
res = 100 ./ length(UT_strain) .* sum( abs(UT(xi, UT_strain) - UT_stress)) ./ max([UT(xi, UT_strain), UT_stress])...
    + 100 ./ length(ET_strain) .* sum( abs(ET(xi, ET_strain) - ET_stress)) ./ max([ET(xi, ET_strain), ET_stress])...
    + 100 ./ length(PS_strain) .* sum( abs(PS(xi, PS_strain) - PS_stress)) ./ max([PS(xi, PS_strain), PS_stress]);

%  normalized mean absolute difference
% res = 100 .* sum( abs(UT(xi, UT_strain) - UT_stress)) ./ length(UT_strain) ./ max( [ sum(abs(UT(xi, UT_strain)))/ length(UT_strain), sum(abs(UT_strain))/length(UT_strain) ] )...
%     + 100 .* sum( abs(ET(xi, ET_strain) - ET_stress)) ./ length(ET_strain) ./ max( [ sum(abs(ET(xi, ET_strain)))/ length(ET_strain), sum(abs(ET_strain))/length(ET_strain) ] )...
%     + 100 .* sum( abs(PS(xi, PS_strain) - PS_stress)) ./ length(PS_strain) ./ max( [ sum(abs(PS(xi, PS_strain)))/ length(PS_strain), sum(abs(PS_strain))/length(PS_strain) ] );

% coefficient of determination
% res = 1.0 - sum((UT(xi, UT_strain) - UT_stress).^2) ./ sum((UT_stress - (max(UT_stress))*ones(size(UT_stress))).^2)...
%     + 1.0 - sum((ET(xi, ET_strain) - ET_stress).^2) ./ sum((ET_stress - (max(ET_stress))*ones(size(ET_stress))).^2)...
%     + 1.0 - sum((PS(xi, PS_strain) - PS_stress).^2) ./ sum((PS_stress - (max(PS_stress))*ones(size(PS_stress))).^2);
% disp(size(res));
end

% Residual
function ff = res(xi, UT_strain, UT_stress, ET_strain, ET_stress, PS_strain, PS_stress, UT, ET, PS)

ff = sum((UT(xi, UT_strain) - UT_stress).^2) ./ length(UT_strain)+ ...
     sum((ET(xi, ET_strain) - ET_stress).^2) ./ length(ET_strain)+ ...
     sum((PS(xi, PS_strain) - PS_stress).^2) ./ length(PS_strain);
end

% Initialize CR Model
function [paras_0, lb, ub, UT, ET, PS] = CR_Model_Init()
lb = [-Inf, -Inf,   0,-Inf, -Inf,   0];
ub = [ Inf,  Inf, Inf, Inf,  Inf, Inf];

paras_0 = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0];

% tool function for generalized strain
term1 = @(x, xdata) 2*x(3)*(xdata.^x(2) - xdata.^(-x(1))) .* ((x(2).*(xdata.^(x(2)-1)) + x(1).*(xdata.^(-x(1)-1)) )  / (x(2)+x(1)).^2);
term2 = @(x, xdata) 2*x(6)*(xdata.^x(5) - xdata.^(-x(4))) .* ((x(5).*(xdata.^(x(5)-1)) + x(4).*(xdata.^(-x(4)-1)) )  / (x(5)+x(4)).^2);

% P_11 of generalized strain
UT = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-1.5)) .* ( term1(x, xdata.^(-0.5)) + term2(x, xdata.^(-0.5)) );
ET = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-3.0)) .* ( term1(x, xdata.^(-2.0)) + term2(x, xdata.^(-2.0)) );
PS = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-2.0)) .* ( term1(x, xdata.^(-1.0)) + term2(x, xdata.^(-1.0)) );
end

% Initialize CR4 Model (4 parameters)
function [paras_0, lb, ub, UT, ET, PS] = CR4_Model_Init()
lb = [-Inf, -Inf, 0, 0];
ub = [Inf, Inf, Inf, Inf];
paras_0 = [1.0, 1.0, 1.0, 1.0];

% tool function for generalized strain
term1 = @(x, xdata) 2*x(3)*(xdata.^x(2) - xdata.^(-x(1))) .* ((x(2).*(xdata.^(x(2)-1)) + x(1).*(xdata.^(-x(1)-1)) )  / (x(2)+x(1)).^2);
term2 = @(x, xdata) 2*x(4).*log(xdata)./xdata;

% P_11 of generalized strain
UT = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-1.5)) .* ( term1(x, xdata.^(-0.5)) + term2(x, xdata.^(-0.5)) );
ET = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-3.0)) .* ( term1(x, xdata.^(-2.0)) + term2(x, xdata.^(-2.0)) );
PS = @(x, xdata) term1(x, xdata) + term2(x, xdata) - (xdata.^(-2.0)) .* ( term1(x, xdata.^(-1.0)) + term2(x, xdata.^(-1.0)) );
end

% Initialize Ogden Model
function [paras_0, lb, ub, UT, ET, PS] = Ogden_Model_Init()
lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf];
ub = [Inf, Inf, Inf, Inf, Inf, Inf];
paras_0 = [1.0, 1.0, 4.0, 2.0, -1.0, -1.0];

UT = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-0.5 * x(2) - 1.0) ) ...
    + x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-0.5 * x(4) - 1.0) )...
    + x(5) * ( xdata .^ (x(6) - 1.0) - xdata .^ (-0.5 * x(6) - 1.0) );

ET = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-2.0 * x(2) - 1.0) ) ...
    + x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-2.0 * x(4) - 1.0) )...
    + x(5) * ( xdata .^ (x(6) - 1.0) - xdata .^ (-2.0 * x(6) - 1.0) );

PS = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-1.0 * x(2) - 1.0) ) ...
    + x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-1.0 * x(4) - 1.0) )...
    + x(5) * ( xdata .^ (x(6) - 1.0) - xdata .^ (-1.0 * x(6) - 1.0) );
end

% Initialize Ogden Model (4 parameters)
function [paras_0, lb, ub, UT, ET, PS] = Ogden4_Model_Init()
lb = [-Inf, -Inf, -Inf, -Inf];
ub = [Inf, Inf, Inf, Inf];
paras_0 = [1.0, 1.0, 1.0, 3.0];

UT = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-0.5 * x(2) - 1.0) ) ...
    + x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-0.5 * x(4) - 1.0) );

ET = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-2.0 * x(2) - 1.0) ) ...
    + x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-2.0 * x(4) - 1.0) );

PS = @(x, xdata) x(1) * ( xdata .^ (x(2) - 1.0) - xdata .^ (-1.0 * x(2) - 1.0) ) ...
    + x(3) * ( xdata .^ (x(4) - 1.0) - xdata .^ (-1.0 * x(4) - 1.0) );
end

% Initialize AB Model
function [paras_0, lb, ub, UT, ET, PS] = AB_Model_Init()
lb = [0, 0];
ub = [Inf, Inf];
paras_0 = [3, 100.0];

UT = @(x, xdata) x(1) .* ( xdata - xdata.^(-2.0) ) .* ( ...
    1 + ( 1.0 ./ ( 5.0 .* x(2) ) ) .* ( xdata.^2 + 2.0.*(xdata.^(-1)) )...
    + ( 33.0 ./ ( 525.0 .* x(2).^2 ) )  .* ( xdata.^2 + 2.0.*(xdata.^(-1)) ).^2 ...
    + ( 76.0 ./ ( 3500.0 .* x(2).^3 ) ) .* ( xdata.^2 + 2.0.*(xdata.^(-1)) ).^3 ...
    + ( 2595.0 ./ ( 336875.0 .* x(2).^4 ) ) .* ( xdata.^2 + 2.0.*(xdata.^(-1)) ).^4 );

ET = @(x, xdata) x(1) .* ( xdata - xdata.^(-5.0) ) .* ( ...
    1 + ( 1.0 ./ ( 5.0 .* x(2) ) ) .* ( 2.*(xdata.^2) + xdata.^(-4) )...
    + ( 33.0 ./ ( 525.0 .* x(2).^2 ) )  .* ( 2.*(xdata.^2) + xdata.^(-4) ).^2 ...
    + ( 76.0 ./ ( 3500.0 .* x(2).^3 ) ) .* ( 2.*(xdata.^2) + xdata.^(-4) ).^3 ...
    + ( 2595.0 ./ ( 336875.0 .* x(2).^4 ) ) .* ( 2.*(xdata.^2) + xdata.^(-4) ).^4 );

PS = @(x, xdata) x(1) .* ( xdata - xdata.^(-3.0) ) .* ( ...
    1 + ( 1.0 ./ ( 5.0 .* x(2) ) ) .* ( xdata.^2 + xdata.^(-2) + 1 )...
    + ( 33.0 ./ ( 525.0 .* x(2).^2 ) )  .* ( xdata.^2 + xdata.^(-2) + 1 ).^2 ...
    + ( 76.0 ./ ( 3500.0 .* x(2).^3 ) ) .* ( xdata.^2 + xdata.^(-2) + 1 ).^3 ...
    + ( 2595.0 ./ ( 336875.0 .* x(2).^4 ) ) .* ( xdata.^2 + xdata.^(-2) + 1 ).^4 );
end

% Initialize MR Model
function [paras_0, lb, ub, UT, ET, PS] = MR_Model_Init()
lb = [0, 0];
ub = [Inf, Inf];
paras_0 = [1, 1];

UT = @(x, xdata) 2.0 .* ( x(1) + x(2)./xdata ) .* ( xdata - xdata.^(-2));
ET = @(x, xdata) 2.0 .* ( x(1) + x(2) .* xdata.^2 ) .* (xdata - xdata.^(-5));
PS = @(x, xdata) 2.0 .* ( x(1) + x(2) ) .* ( xdata - xdata.^(-3));
end
