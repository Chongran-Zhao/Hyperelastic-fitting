function plot_results_James(paras, ...
                      F_UT_list, P_UT_list, ...
                      F_BT_list_1, relative_P1_BT_1, ...
                      F_BT_list_2, relative_P1_BT_2, ...
                      F_BT_list_3, relative_P1_BT_3, ...
                      F_BT_list_4, relative_P1_BT_4, ...
                      F_BT_list_5, relative_P1_BT_5, ...
                      F_BT_list_6, relative_P1_BT_6, ...
                      F_BT_list_7, relative_P1_BT_7, ...
                      F_BT_list_8, relative_P1_BT_8, ...
                      F_BT_list_9, relative_P1_BT_9, ...
                      F_BT_list_10, relative_P1_BT_10, ...
                      F_BT_list_11, relative_P1_BT_11, ...
                      F_BT_list_12, relative_P1_BT_12, ...
                      F_BT_list_13, relative_P1_BT_13, ...
                      F_BT_list_14, relative_P1_BT_14)

P_UT_fit = get_P_ij_list(1, 1, paras, F_UT_list);
UT_lambda_1 = zeros(length(P_UT_list), 1);
UT_lambda_1(:) = F_UT_list(1,1,:);

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

relative_P1_BT_fit_5 = get_P_ij_list(1, 1, paras, F_BT_list_5);
BT_lambda_5 = zeros(length(relative_P1_BT_fit_5), 1);
BT_lambda_5(:) = F_BT_list_5(1,1,:);

relative_P1_BT_fit_6 = get_P_ij_list(1, 1, paras, F_BT_list_6);
BT_lambda_6 = zeros(length(relative_P1_BT_fit_6), 1);
BT_lambda_6(:) = F_BT_list_6(1,1,:);

relative_P1_BT_fit_7 = get_P_ij_list(1, 1, paras, F_BT_list_7);
BT_lambda_7 = zeros(length(relative_P1_BT_fit_7), 1);
BT_lambda_7(:) = F_BT_list_7(1,1,:);

relative_P1_BT_fit_8 = get_P_ij_list(2, 2, paras, F_BT_list_8);
BT_lambda_8 = zeros(length(relative_P1_BT_fit_8), 1);
BT_lambda_8(:) = F_BT_list_8(1,1,:);

relative_P1_BT_fit_9 = get_P_ij_list(2, 2, paras, F_BT_list_9);
BT_lambda_9 = zeros(length(relative_P1_BT_fit_9), 1);
BT_lambda_9(:) = F_BT_list_9(1,1,:);

relative_P1_BT_fit_10 = get_P_ij_list(2, 2, paras, F_BT_list_10);
BT_lambda_10 = zeros(length(relative_P1_BT_fit_10), 1);
BT_lambda_10(:) = F_BT_list_10(1,1,:);

relative_P1_BT_fit_11 = get_P_ij_list(2, 2, paras, F_BT_list_11);
BT_lambda_11 = zeros(length(relative_P1_BT_fit_11), 1);
BT_lambda_11(:) = F_BT_list_11(1,1,:);

relative_P1_BT_fit_12 = get_P_ij_list(2, 2, paras, F_BT_list_12);
BT_lambda_12 = zeros(length(relative_P1_BT_fit_12), 1);
BT_lambda_12(:) = F_BT_list_12(1,1,:);

relative_P1_BT_fit_13 = get_P_ij_list(2, 2, paras, F_BT_list_13);
BT_lambda_13 = zeros(length(relative_P1_BT_fit_13), 1);
BT_lambda_13(:) = F_BT_list_13(1,1,:);

relative_P1_BT_fit_14 = get_P_ij_list(2, 2, paras, F_BT_list_14);
BT_lambda_14 = zeros(length(relative_P1_BT_fit_14), 1);
BT_lambda_14(:) = F_BT_list_14(1,1,:);

figure;

hold on;
plot(UT_lambda_1, P_UT_fit, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(UT_lambda_1, P_UT_list, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
hold off;

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
plot(BT_lambda_3, relative_P1_BT_3, 'go', 'MarkerFaceColor', [0.0, 0.6, 0.0], 'MarkerSize', 12); 
plot(BT_lambda_4, relative_P1_BT_4, 'bo', 'MarkerFaceColor', [0.00, 0.45, 0.74], 'MarkerSize', 12); 
plot(BT_lambda_5, relative_P1_BT_5, 'mo', 'MarkerFaceColor', [0.5, 0.0, 0.5], 'MarkerSize', 12); 
plot(BT_lambda_6, relative_P1_BT_6, 'co', 'MarkerFaceColor', [0.0, 0.75, 0.75], 'MarkerSize', 12); 
plot(BT_lambda_7, relative_P1_BT_7, 'yo', 'MarkerFaceColor', [1.0, 0.75, 0.0], 'MarkerSize', 12); 

plot(BT_lambda_1, relative_P1_BT_fit_1, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_2, relative_P1_BT_fit_2, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_3, relative_P1_BT_fit_3, 'Color', [0.0, 0.6, 0.0], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_4, relative_P1_BT_fit_4, 'Color', [0.00, 0.45, 0.74], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_5, relative_P1_BT_fit_5, 'Color', [0.5, 0.0, 0.5], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_6, relative_P1_BT_fit_6, 'Color', [0.0, 0.75, 0.75], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_7, relative_P1_BT_fit_7, 'Color', [1.0, 0.75, 0.0], 'LineWidth', 2.5, 'LineStyle', '-');

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

l = legend( '$\lambda_2 = 3.5$', '$\lambda_2 = 3.0$', '$\lambda_2 = 2.5$', ...
            '$\lambda_2 = 2.0$', '$\lambda = 1.7$', '$\lambda=1.5$', '$\lambda=1.3$' );
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

plot(BT_lambda_8, relative_P1_BT_8, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(BT_lambda_9, relative_P1_BT_9, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(BT_lambda_10, relative_P1_BT_10, 'go', 'MarkerFaceColor', [0.0, 0.6, 0.0], 'MarkerSize', 12); 
plot(BT_lambda_11, relative_P1_BT_11, 'bo', 'MarkerFaceColor', [0.00, 0.45, 0.74], 'MarkerSize', 12); 
plot(BT_lambda_12, relative_P1_BT_12, 'mo', 'MarkerFaceColor', [0.5, 0.0, 0.5], 'MarkerSize', 12); 
plot(BT_lambda_13, relative_P1_BT_13, 'co', 'MarkerFaceColor', [0.0, 0.75, 0.75], 'MarkerSize', 12); 
plot(BT_lambda_14, relative_P1_BT_14, 'yo', 'MarkerFaceColor', [1.0, 0.75, 0.0], 'MarkerSize', 12); 

plot(BT_lambda_8, relative_P1_BT_fit_8, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_9, relative_P1_BT_fit_9, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_10, relative_P1_BT_fit_10, 'Color', [0.0, 0.6, 0.0], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_11, relative_P1_BT_fit_11, 'Color', [0.00, 0.45, 0.74], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_12, relative_P1_BT_fit_12, 'Color', [0.5, 0.0, 0.5], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_13, relative_P1_BT_fit_13, 'Color', [0.0, 0.75, 0.75], 'LineWidth', 2.5, 'LineStyle', '-');
plot(BT_lambda_14, relative_P1_BT_fit_14, 'Color', [1.0, 0.75, 0.0], 'LineWidth', 2.5, 'LineStyle', '-');

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

l = legend( '$\lambda_2 = 3.5$', '$\lambda_2 = 3.0$', '$\lambda_2 = 2.5$', ...
            '$\lambda_2 = 2.0$', '$\lambda = 1.7$', '$\lambda=1.5$', '$\lambda=1.3$' );
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