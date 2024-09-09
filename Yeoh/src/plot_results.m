function plot_results(xi, F_list_1, P_exp_1,...
                          F_list_2, P_exp_2,...
                          F_list_3, P_exp_3)

P_fit_1 = get_P_ij_list(1, 1, xi, F_list_1);
P_fit_2 = get_P_ij_list(1, 1, xi, F_list_2);
P_fit_3 = get_P_ij_list(1, 1, xi, F_list_3);

lambda_1 = zeros(length(P_exp_1), 1);
lambda_2 = zeros(length(P_exp_2), 1);
lambda_3 = zeros(length(P_exp_3), 1);

lambda_1(:) = F_list_1(1,1,:);
lambda_2(:) = F_list_2(1,1,:);
lambda_3(:) = F_list_3(1,1,:);

figure;
ax = axes('Position', [0.1 0.4 0.8 0.5], 'Box', 'on');
plot(ax, lambda_1, P_exp_1, 'Color', '#003f5c', 'Marker', 'o', 'MarkerFaceColor', '#003f5c', 'MarkerSize', 12, 'LineStyle', 'none');
hold(ax, 'on');
plot(ax, lambda_1, P_fit_1, 'linewidth', 3.0, 'Color', '#003f5c', 'LineStyle', '-');
hold(ax, 'on');
plot(ax, lambda_2, P_exp_2, 'Color', '#58508d', 'Marker', 'o', 'MarkerFaceColor', '#58508d', 'MarkerSize', 12, 'LineStyle', 'none');
hold(ax, 'on');
plot(ax, lambda_2, P_fit_2, 'linewidth', 3.0, 'Color', '#58508d', 'LineStyle', '-');
hold(ax, 'on');
plot(ax, lambda_3, P_exp_3, 'Color', '#bc5090', 'Marker', 'o', 'MarkerFaceColor', '#bc5090', 'MarkerSize', 12, 'LineStyle', 'none');
hold(ax, 'on');
plot(ax, lambda_3, P_fit_3, 'linewidth', 3.0, 'Color', '#bc5090', 'LineStyle', '-');
xlabel(ax, 'Stretch', 'interpreter', 'latex', 'FontSize', 30, 'FontWeight', 'bold', 'FontName', 'Helvetica');
ylabel(ax, 'Nominal stress', 'interpreter', 'latex', 'FontSize', 30, 'FontWeight', 'bold', 'FontName', 'Helvetica');
set(ax, 'TickDir', 'out', ...
    'TickLength', [.02 .02], ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'YGrid', 'on', ...
    'XGrid', 'on', ...
    'XColor', [0 0 0], ...
    'YColor', [0 0 0], ...
    'LineWidth', 2, ...
    'FontSize', 25, 'FontWeight', 'bold');
l = legend(ax, 'exp-UT', 'pre-UT',...
    'exp-ET', 'pre-ET',...
    'exp-PS', 'pre-PS',...
    'location', 'northwest', 'Orientation', 'horizontal');
set(l, 'interpreter', 'latex', 'fontsize', 25, 'box', 'off', 'FontWeight', 'bold', 'FontName', 'Helvetica', 'NumColumns', 4);

X = 40.0;
Y = 40.0;
xMargin = 3;
yMargin = 3;
xSize = X - 2 * xMargin;
ySize = Y - 2 * yMargin;
set(gcf, 'Units', 'centimeters', 'Position', [5 5 xSize ySize]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [X Y]);
set(gcf, 'PaperPosition', [xMargin yMargin xSize ySize]);
set(gcf, 'PaperOrientation', 'portrait');

x_location = 7.0;
y_location = 0.1;
delta_y = 0.5;
% print R^2
R_square = (get_R_square(P_exp_1, P_fit_1) + get_R_square(P_exp_2, P_fit_2) + get_R_square(P_exp_3, P_fit_3)) / 3.0;
text_R_square = sprintf('$R^2=%.4g$ \\%%', 100*R_square);
text(x_location, y_location, text_R_square, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', ...
    'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');

% print NMAD
y_location = y_location + delta_y;
NMAD = (get_NMAD(P_fit_1, P_exp_1) + get_NMAD(P_fit_2, P_exp_2) + get_NMAD(P_fit_3, P_exp_3)) / 3.0;
text_NMAD = sprintf('$\\mathrm{NMAD}=%.4g$ \\%%', NMAD);
text(x_location, y_location, text_NMAD, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', ...
    'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');

% print MSD
y_location = y_location + delta_y;
MSD = (get_MSD(P_fit_1, P_exp_1) + get_MSD(P_fit_1, P_exp_1) + get_MSD(P_fit_1, P_exp_1)) / 3.0;
text_MSD = sprintf('$\\mathrm{MSD} = %.4g$', MSD);
text(x_location, y_location, text_MSD, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', ...
    'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');

% print parameters
text_xi = cell(length(xi), 1);
x_location = 2.0;
y_location = -2.0;
delta_x = 2.0;
for ii = 1:length(xi)
    text_xi{ii} = sprintf('$C_%d = %.4g$', ii, xi(ii));
    text(x_location, y_location, text_xi{ii}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Interpreter', 'latex', ...
        'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Helvetica');
    x_location = x_location + delta_x;
end
end