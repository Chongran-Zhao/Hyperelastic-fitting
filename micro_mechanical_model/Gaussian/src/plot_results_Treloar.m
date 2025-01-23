function plot_results_Treloar(paras, ...
                      F_UT_list, P_UT_list, ...
                      F_ET_list, P_ET_list, ...
                      F_PS_list, P_PS_list)

P_UT_fit = get_P_ij_list(1, 1, paras, F_UT_list);
UT_lambda_1 = zeros(length(P_UT_list), 1);
UT_lambda_1(:) = F_UT_list(1,1,:);

P_ET_fit = get_P_ij_list(1, 1, paras, F_ET_list);
ET_lambda_1 = zeros(length(P_ET_list), 1);
ET_lambda_1(:) = F_ET_list(1,1,:);

P_PS_fit = get_P_ij_list(1, 1, paras, F_PS_list);
PS_lambda_1 = zeros(length(P_PS_list), 1);
PS_lambda_1(:) = F_PS_list(1,1,:);

figure;

hold on;
plot(UT_lambda_1, P_UT_list, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 
plot(UT_lambda_1, P_UT_fit, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
plot(ET_lambda_1, P_ET_list, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(ET_lambda_1, P_ET_fit, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(PS_lambda_1, P_PS_list, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 
plot(PS_lambda_1, P_PS_fit, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');
hold off;
l = legend('UT-fitting', 'ET-prediction', 'PS-prediction' );
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

end