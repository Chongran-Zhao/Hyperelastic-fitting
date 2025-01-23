function plot_results_Kawamura(paras, ...
                      F_UT_list, P_UT_list, ...
                      F_UE_list, P_UE_list, ...
                      F_ET_list, P_ET_list, ...
                      F_BT_list_1, relative_P1_BT_1, ...
                      F_BT_list_2, relative_P1_BT_2, ...
                      F_BT_list_3, relative_P1_BT_3, ...
                      F_BT_list_4, relative_P1_BT_4, ...
                      F_BT_list_5, relative_P1_BT_5, ...
                      F_BT_list_6, relative_P1_BT_6, ...
                      F_BT_list_7, relative_P1_BT_7, ...
                      F_BT_list_8, relative_P1_BT_8)

P_UT_fit = get_P_ij_list(1, 1, paras, F_UT_list);
UT_lambda_1 = zeros(length(P_UT_list), 1);
UT_lambda_1(:) = F_UT_list(1,1,:);

P_UE_fit = get_P_ij_list(1, 1, paras, F_UE_list);
UE_lambda_1 = zeros(length(P_UE_list), 1);
UE_lambda_1(:) = F_UE_list(1,1,:);

P_ET_fit = get_P_ij_list(1, 1, paras, F_ET_list);
ET_lambda_1 = zeros(length(P_ET_list), 1);
ET_lambda_1(:) = F_ET_list(1,1,:);

relative_P1_BT_fit_1 = get_P_ij_list(1, 1, paras, F_BT_list_1);
BT_lambda_1 = zeros(length(relative_P1_BT_fit_1), 1);
BT_lambda_1(:) = F_BT_list_1(1,1,:);

relative_P1_BT_fit_2 = get_P_ij_list(1, 1, paras, F_BT_list_2);
BT_lambda_2 = zeros(length(relative_P1_BT_fit_2), 1);
BT_lambda_2(:) = F_BT_list_2(1,1,:);

relative_P1_BT_fit_3 = get_P_ij_list(1, 1, paras, F_BT_list_3);
BT_lambda_3 = zeros(length(relative_P1_BT_fit_3), 1);
BT_lambda_3(:) = F_BT_list_3(1,1,:);

relative_P1_BT_fit_4 = get_P_ij_list(1, 1, paras, F_BT_list_4);
BT_lambda_4 = zeros(length(relative_P1_BT_fit_4), 1);
BT_lambda_4(:) = F_BT_list_4(1,1,:);

relative_P1_BT_fit_5 = get_P_ij_list(2, 2, paras, F_BT_list_5);
BT_lambda_5 = zeros(length(relative_P1_BT_fit_5), 1);
BT_lambda_5(:) = F_BT_list_5(1,1,:);

relative_P1_BT_fit_6 = get_P_ij_list(2, 2, paras, F_BT_list_6);
BT_lambda_6 = zeros(length(relative_P1_BT_fit_6), 1);
BT_lambda_6(:) = F_BT_list_6(1,1,:);

relative_P1_BT_fit_7 = get_P_ij_list(2, 2, paras, F_BT_list_7);
BT_lambda_7 = zeros(length(relative_P1_BT_fit_7), 1);
BT_lambda_7(:) = F_BT_list_7(1,1,:);

relative_P1_BT_fit_8 = get_P_ij_list(2, 2, paras, F_BT_list_8);
BT_lambda_8 = zeros(length(relative_P1_BT_fit_8), 1);
BT_lambda_8(:) = F_BT_list_8(1,1,:);

figure;

hold on;
plot(UT_lambda_1, P_UT_fit, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(UE_lambda_1, P_UE_fit, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(ET_lambda_1, P_ET_fit, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');
plot(UT_lambda_1, P_UT_list, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(UE_lambda_1, P_UE_list, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(ET_lambda_1, P_ET_list, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 
hold off;
l = legend('UT-fitting', 'UC-prediction', 'ET-prediction' );
set( l, 'interpreter','latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', 'Orientation', ...
     'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 3 );
l.ItemTokenSize=[45,45];

title('Fitting of Uniaxial Tensile Test', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
hXLabel = xlabel('$\lambda_1$', 'interpreter', 'latex');
hYLabel = ylabel('$P_1$ (MPa)', 'interpreter', 'latex');
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

figure;
title('Prediction of Biaxial Tensile Test', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
hold on;
plot(BT_lambda_1, relative_P1_BT_1, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(BT_lambda_2, relative_P1_BT_2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(BT_lambda_3, relative_P1_BT_3, 'o', 'MarkerFaceColor', '0.0, 0.6, 0.0', 'MarkerSize', 12); 
plot(BT_lambda_4, relative_P1_BT_4, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 

plot(BT_lambda_1, relative_P1_BT_fit_1, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_2, relative_P1_BT_fit_2, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_3, relative_P1_BT_fit_3, 'Color', '0.0, 0.6, 0.0', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_4, relative_P1_BT_fit_4, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');
hold off;

hXLabel = xlabel('$\lambda_1$', 'interpreter', 'latex');
hYLabel = ylabel('$P_1$ (MPa)', 'interpreter', 'latex');
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

l = legend( '$\lambda_2 = 1.7$', '$\lambda_2 = 1.5$', '$\lambda_2 = 1.3$', ...
            '$\lambda_2 = 1.1$' );
set( l, 'interpreter','latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', 'Orientation', ...
     'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 4 );
l.ItemTokenSize=[45,45];

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

figure;
title('Prediction of Biaxial Tensile Test', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
hold on;
plot(BT_lambda_5, relative_P1_BT_5, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(BT_lambda_6, relative_P1_BT_6, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(BT_lambda_7, relative_P1_BT_7, 'o', 'MarkerFaceColor', '0.0, 0.6, 0.0', 'MarkerSize', 12); 
plot(BT_lambda_8, relative_P1_BT_8, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 

plot(BT_lambda_5, relative_P1_BT_fit_5, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_6, relative_P1_BT_fit_6, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_7, relative_P1_BT_fit_7, 'Color', '0.0, 0.6, 0.0', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_8, relative_P1_BT_fit_8, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');
hold off;

hXLabel = xlabel('$\lambda_1$', 'interpreter', 'latex');
hYLabel = ylabel('$P_2$ (MPa)', 'interpreter', 'latex');
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

l = legend( '$\lambda_2 = 1.7$', '$\lambda_2 = 1.5$', '$\lambda_2 = 1.3$', ...
            '$\lambda_2 = 1.1$' );
set( l, 'interpreter','latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', 'Orientation', ...
     'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 4 );
l.ItemTokenSize=[45,45];

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