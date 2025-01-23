function plot_results_Katashima(paras,...
             F_ET_list, P_ET_list, ...
             F_BT_list_1, relative_P1_BT_1, ...
             F_BT_list_2, relative_P1_BT_2, ...
             F_PS_list_1, relative_P1_PS_1, ...
             F_PS_list_2, relative_P1_PS_2);

P_ET_fit = get_P_ij_list(1, 1, paras, F_ET_list);
ET_lambda_1 = zeros(length(P_ET_list), 1);
ET_lambda_1(:) = F_ET_list(1,1,:);

relative_P1_BT_fit_1 = get_P_ij_list(1, 1, paras, F_BT_list_1);
BT_lambda_1 = zeros(length(relative_P1_BT_fit_1), 1);
BT_lambda_1(:) = F_BT_list_1(1,1,:);

relative_P1_BT_fit_2 = get_P_ij_list(2, 2, paras, F_BT_list_2);
BT_lambda_2 = zeros(length(relative_P1_BT_fit_2), 1);
BT_lambda_2(:) = F_BT_list_2(1,1,:);

relative_P1_PS_fit_1 = get_P_ij_list(1, 1, paras, F_PS_list_1);
PS_lambda_1 = zeros(length(relative_P1_PS_fit_1), 1);
PS_lambda_1(:) = F_PS_list_1(1,1,:);

relative_P1_PS_fit_2 = get_P_ij_list(2, 2, paras, F_PS_list_2);
PS_lambda_2 = zeros(length(relative_P1_PS_fit_2), 1);
PS_lambda_2(:) = F_PS_list_2(1,1,:);

figure;

hold on;
plot(ET_lambda_1, P_ET_fit, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');

plot(BT_lambda_1, relative_P1_BT_1, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(BT_lambda_2, relative_P1_BT_2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(PS_lambda_1, relative_P1_PS_1, 'o', 'MarkerFaceColor', '0.0, 0.6, 0.0', 'MarkerSize', 12); 
plot(PS_lambda_2, relative_P1_PS_2, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 

plot(ET_lambda_1, P_ET_list, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 

plot(BT_lambda_1, relative_P1_BT_fit_1, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_2, relative_P1_BT_fit_2, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(PS_lambda_1, relative_P1_PS_fit_1, 'Color', '0.0, 0.6, 0.0', 'LineWidth', 2.5, 'LineStyle', '-');
plot(PS_lambda_2, relative_P1_PS_fit_2, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');
hold off;
l = legend('ET-fitting', 'BT-prediction-$P_1$', 'BT-prediction-$P_2$', 'PS-prediction-$P_1$', 'PS-prediction-$P_2$' );
set( l, 'interpreter','latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', 'Orientation', ...
     'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 3 );
l.ItemTokenSize=[45,45];

title('Fitting of Equi-Biaxial Tensile Test', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
hXLabel = xlabel('$\lambda_1$', 'interpreter', 'latex');
hYLabel = ylabel('$P_1$ (kPa)', 'interpreter', 'latex');
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
X = 42.0;
Y = X * 0.9;
xMargin = 3;
yMargin = 3;
xSize = X - 2 * xMargin;
ySize = Y - 2 * yMargin;
set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[X Y]);
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize]);
set(gcf, 'PaperOrientation','portrait');

% print NMAD
% NMAD = (get_NMAD(P_UT_fit, P_UT_list) + get_NMAD(P_ET_pre, P_ET_list) + get_NMAD(P_PS_pre, P_PS_list)) / 3.0;
% text_NMAD = sprintf('$\\mathrm{NMAD}=%.4g \\%%$', NMAD);
% text(max([max(lambda_1),max(lambda_2),max(lambda_3)])-0.5, -1.5, text_NMAD, ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 50, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');


end