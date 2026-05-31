function figureHandle = plot_simultaneous_fit(models, cases)
% Plot simultaneous fitting results in one shared coordinate system.

caseList = normalize_cases(cases);
figureHandle = figure('Name', 'Simultaneous fitting results', 'Color', 'w');
axesHandle = axes(figureHandle);
hold(axesHandle, 'on');

colors = scientific_colors();
markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
seriesIndex = 0;
experimentHandles = [];
fitHandles = [];
experimentLabels = {};
fitLabels = {};
evaluationSeries = empty_series();

for caseIndex = 1:length(caseList)
    caseData = caseList{caseIndex};
    seriesList = collect_case_series(models, caseData);

    for ii = 1:length(seriesList)
        seriesIndex = seriesIndex + 1;
        color = colors(wrap_index(seriesIndex, size(colors, 1)), :);
        marker = markers{wrap_index(seriesIndex, length(markers))};
        series = seriesList(ii);

        experimentHandles(end + 1) = plot_experiment_series(series, color, marker);
        fitHandles(end + 1) = plot_model_series(series, color);
        experimentLabels{end + 1} = sprintf('%s exp', series.label);
        fitLabels{end + 1} = sprintf('%s fit', series.label);
        evaluationSeries(end + 1) = series;
    end
end

set_axes_style(axesHandle);
xlabel('$\lambda_1$', 'Interpreter', 'latex', ...
    'FontName', 'Times New Roman');
ylabel('$P_{11}$ (MPa)', 'Interpreter', 'latex', ...
    'FontName', 'Times New Roman');
Add_evaluation_text(axesHandle, Evaluation_metrics(evaluationSeries));
legendHandle = legend([experimentHandles, fitHandles], ...
    [experimentLabels, fitLabels], ...
    'Location', 'southoutside', ...
    'Orientation', 'horizontal', ...
    'NumColumns', max(1, seriesIndex), ...
    'Box', 'off', ...
    'FontName', 'Times New Roman');
legendHandle.ItemTokenSize = [22, 12];
end

function caseList = normalize_cases(cases)
if iscell(cases)
    caseList = cases;
elseif isstruct(cases)
    caseList = num2cell(cases);
else
    error('plot_simultaneous_fit:InvalidCases', ...
        'cases must be a cell array or struct array.');
end
end

function seriesList = collect_case_series(models, caseData)
if strcmpi(caseData.stress_type, 'scalar')
    seriesList = scalar_case_series(models, caseData);
elseif strcmpi(caseData.stress_type, 'relative_cauchy')
    seriesList = relative_cauchy_case_series(models, caseData);
else
    seriesList = tensor_case_series(models, caseData);
end
end

function seriesList = scalar_case_series(models, caseData)
[row, col] = component_indices(caseData.stress_component);

seriesList = empty_series();
seriesList(1).lambda = caseData.lambda(:, 1);
seriesList(1).exp = caseData.P_exp(:);
seriesList(1).fit = fitted_PK1_component(models, caseData, row, col);
seriesList(1).label = sprintf('%s %s', caseData.name, caseData.stress_component);
end

function seriesList = tensor_case_series(models, caseData)
components = {'P11', 'P22', 'P33', 'P12', 'P21'};
seriesList = empty_series();

for componentIndex = 1:length(components)
    [row, col] = component_indices(components{componentIndex});
    mask = reshape(caseData.P_mask(row, col, :), [], 1);
    if ~any(mask)
        continue;
    end

    P_exp_all = reshape(caseData.P_list(row, col, :), [], 1);
    P_fit_all = fitted_PK1_component(models, caseData, row, col);

    series.lambda = caseData.lambda(mask, 1);
    series.exp = P_exp_all(mask);
    series.fit = P_fit_all(mask);
    series.label = sprintf('%s %s', caseData.name, components{componentIndex});
    seriesList(end + 1) = series;
end
end

function seriesList = relative_cauchy_case_series(models, caseData)
seriesList = empty_series();
seriesList(1).lambda = caseData.lambda(:, 1);
seriesList(1).exp = caseData.sigma_exp(:);
seriesList(1).fit = fitted_relative_cauchy(models, caseData);
seriesList(1).label = sprintf('%s sigma11-sigma22', caseData.name);
end

function seriesList = empty_series()
seriesList = struct('lambda', {}, 'exp', {}, 'fit', {}, 'label', {});
end

function plotHandle = plot_experiment_series(series, color, marker)
plotHandle = plot(series.lambda, series.exp, marker, ...
    'LineStyle', 'none', ...
    'MarkerSize', 5.5, ...
    'LineWidth', 1.1, ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerEdgeColor', color, ...
    'DisplayName', sprintf('%s exp', series.label));
end

function plotHandle = plot_model_series(series, color)
[lambdaSorted, fitSorted] = sort_curve(series.lambda, series.fit);
plotHandle = plot(lambdaSorted, fitSorted, '-', ...
    'Color', color, ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('%s fit', series.label));
end

function values = fitted_PK1_component(models, caseData, row, col)
values = zeros(caseData.num_points, 1);

for pointIndex = 1:caseData.num_points
    F = caseData.F_list(:, :, pointIndex);
    P = models.P(F);
    values(pointIndex) = P(row, col);
end
end

function values = fitted_relative_cauchy(models, caseData)
values = zeros(caseData.num_points, 1);

for pointIndex = 1:caseData.num_points
    F = caseData.F_list(:, :, pointIndex);
    P = models.P(F);
    values(pointIndex) = relative_cauchy_stress(P, F, ...
        caseData.relative_cauchy_component);
end
end

function value = relative_cauchy_stress(P, F, name)
J = det(F);
sigma = P * F' ./ J;

switch upper(char(name))
    case 'SIGMA11_MINUS_SIGMA22'
        value = sigma(1, 1) - sigma(2, 2);
    otherwise
        error('plot_simultaneous_fit:UnsupportedRelativeCauchyComponent', ...
            'Unsupported relative Cauchy component: %s.', name);
end
end

function colors = scientific_colors()
colors = [ ...
    0.00, 0.45, 0.74; ...
    0.85, 0.33, 0.10; ...
    0.47, 0.67, 0.19; ...
    0.49, 0.18, 0.56; ...
    0.30, 0.75, 0.93; ...
    0.93, 0.69, 0.13; ...
    0.64, 0.08, 0.18; ...
    0.20, 0.20, 0.20; ...
    0.10, 0.55, 0.40];
end

function index = wrap_index(index, count)
index = mod(index - 1, count) + 1;
end

function [lambdaSorted, valueSorted] = sort_curve(lambda, value)
[lambdaSorted, order] = sort(lambda);
valueSorted = value(order);
end

function set_axes_style(axesHandle)
set(axesHandle, ...
    'Box', 'off', ...
    'LineWidth', 1.1, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'TickDir', 'out', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Layer', 'top');
grid(axesHandle, 'on');
axesHandle.GridAlpha = 0.14;
axesHandle.MinorGridAlpha = 0.08;
end

function [row, col] = component_indices(name)
switch upper(char(name))
    case 'P11'
        row = 1; col = 1;
    case 'P22'
        row = 2; col = 2;
    case 'P33'
        row = 3; col = 3;
    case 'P12'
        row = 1; col = 2;
    case 'P21'
        row = 2; col = 1;
    otherwise
        error('plot_simultaneous_fit:UnsupportedComponent', ...
            'Unsupported PK1 component: %s.', name);
end
end
