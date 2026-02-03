function plot_results(names, paras,...
    F_list_1, P_exp_1,...
    F_list_2, P_exp_2,...
    F_list_3, P_exp_3)
mu = [];
m = [];
n = [];
for ii = 1:length(names)
    mu = [mu, paras(3*ii-2)];
    m = [m, paras(3*ii-1)];
    n = [n, paras(3*ii)];
end
P_fit_1 = get_P_ij_list(1, 1, names, paras, F_list_1);
P_fit_2 = get_P_ij_list(1, 1, names, paras, F_list_2);
P_fit_3 = get_P_ij_list(1, 1, names, paras, F_list_3);

lambda_1 = zeros(length(P_exp_1), 1);
lambda_2 = zeros(length(P_exp_2), 1);
lambda_3 = zeros(length(P_exp_3), 1);

lambda_1(:) = F_list_1(1,1,:);
lambda_2(:) = F_list_2(1,1,:);
lambda_3(:) = F_list_3(1,1,:);

figure;
hold on;
plot(lambda_1, P_exp_1, 'o', 'MarkerFaceColor', '0.00,0.45,0.74', 'MarkerSize', 12); 
plot(lambda_2, P_exp_2, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12); 
plot(lambda_3, P_exp_3, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12); 

plot(lambda_1, P_fit_1, 'Color', '0.00,0.45,0.74', 'LineWidth', 2.5, 'LineStyle', '-');
plot(lambda_2, P_fit_2, 'Color', 'r', 'LineWidth', 2.5, 'LineStyle', '-');
plot(lambda_3, P_fit_3, 'Color', 'k', 'LineWidth', 2.5, 'LineStyle', '-');
hold off;
% ylim( [0, 6] );
% xlim( [1, 8] );
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

% l = legend( 'UT Experimental', 'ET Experimental', 'PS Experimental', ...
%             'UT Simulated', 'ET Simulated', 'PS Simulated' );
% set( l, 'interpreter','latex', 'fontsize', 25, 'box', 'off', 'location', 'SouthOutside', 'Orientation', ...
%      'horizontal', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 3 );
% l.ItemTokenSize=[45,45];

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

% x_location = 7.0;
% y_location = 0.1;
% delta_y = 0.5;
% print R^2
% R_square = (get_R_square(P_exp_1, P_fit_1) + get_R_square(P_exp_2, P_fit_2) + get_R_square(P_exp_3, P_fit_3)) / 3.0;
% text_R_square = sprintf('$R^2=%.4g \\%%$', 100*R_square);
% text(x_location, y_location, text_R_square, ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');

% print NMAD
NMAD = (get_NMAD(P_fit_1, P_exp_1) + get_NMAD(P_fit_2, P_exp_2) + get_NMAD(P_fit_3, P_exp_3)) / 3.0;
text_NMAD = sprintf('$\\mathrm{NMAD}=%.4g \\%%$', NMAD);
text(max([max(lambda_1),max(lambda_2),max(lambda_3)])-1.3, 0.0, text_NMAD, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', ...
    'FontSize', 50, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');

% print MSD
% y_location = y_location + delta_y;
% MSD = (get_MSD(P_fit_1, P_exp_1) + get_MSD(P_fit_2, P_exp_2) + get_MSD(P_fit_3, P_exp_3)) / 3.0;
% text_MSD = sprintf('$\\mathrm{MSD} = %.4g$', MSD);
% text(x_location, y_location, text_MSD, ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'bottom', ...
%     'Interpreter', 'latex', ...
%     'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');

% print parameters
% text_mu = cell(length(mu), 1);
% text_m = cell(length(m), 1);
% text_n = cell(length(n), 1);
% y_location = -2.0;
% delta_x = 2.0;
% delta_y = -1.0;
% for ii = 1:length(mu)
%     x_location = 2.0;
%     text_mu{ii} = sprintf('$\\mu_%d = %.4g$', ii, mu(ii));
%     text(x_location, y_location, text_mu{ii}, ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'Interpreter', 'latex', ...
%         'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');
%     x_location = x_location + delta_x;
% 
%     text_m{ii} = sprintf('$m_%d = %.4g$', ii, m(ii));
%     text(x_location, y_location, text_m{ii}, ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'Interpreter', 'latex', ...
%         'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');
%     x_location = x_location + delta_x;
% 
%     text_n{ii} = sprintf('$n_%d = %.4g$', ii, n(ii));
%     text(x_location, y_location, text_n{ii}, ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'Interpreter', 'latex', ...
%         'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');
%     y_location = y_location + delta_y;
% end
end