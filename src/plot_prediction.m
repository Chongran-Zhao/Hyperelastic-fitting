function figureHandle = plot_prediction(models, predictionCases)
% Plot prediction results using fitted material parameters.

caseList = normalize_cases(predictionCases);
if isempty(caseList)
    figureHandle = [];
    return;
end

panels = build_prediction_panels(models, caseList);
if isempty(panels)
    figureHandle = [];
    return;
end

numPanels = length(panels);

figureHandle = figure('Name', 'Prediction results', 'Color', 'w');
[numRows, numCols] = prediction_layout_size(numPanels);
layoutHandle = tiledlayout(numRows, numCols, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
sharedExperimentHandles = [];
sharedPredictionHandles = [];
sharedExperimentLabels = {};
sharedPredictionLabels = {};
sharedLegendKeys = {};
legendAxesHandle = [];

for panelIndex = 1:numPanels
    axesHandle = nexttile(layoutHandle);
    if isempty(legendAxesHandle)
        legendAxesHandle = axesHandle;
    end
    hold(axesHandle, 'on');

    [experimentHandles, predictionHandles, experimentLabels, predictionLabels, legendKeys] = ...
        plot_panel(panels(panelIndex), axesHandle);
    set_axes_style(axesHandle);

    xlabel('$\lambda_1$', 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman');
    ylabel(panels(panelIndex).stress_label, 'Interpreter', 'latex', ...
        'FontName', 'Times New Roman');
    add_evaluation_text(axesHandle, ...
        evaluation_metrics(panels(panelIndex).seriesList));

    [sharedExperimentHandles, sharedPredictionHandles, ...
        sharedExperimentLabels, sharedPredictionLabels, sharedLegendKeys] = ...
        append_shared_legend_entries(sharedExperimentHandles, ...
        sharedPredictionHandles, sharedExperimentLabels, ...
        sharedPredictionLabels, sharedLegendKeys, experimentHandles, ...
        predictionHandles, experimentLabels, predictionLabels, legendKeys);
end

if ~isempty(sharedExperimentHandles)
    legendHandle = legend(legendAxesHandle, ...
        [sharedExperimentHandles, sharedPredictionHandles], ...
        [sharedExperimentLabels, sharedPredictionLabels], ...
        'Orientation', 'horizontal', ...
        'NumColumns', length(sharedExperimentHandles), ...
        'Box', 'off', ...
        'Interpreter', 'tex', ...
        'FontName', 'Times New Roman');
    legendHandle.ItemTokenSize = [22, 12];
    legendHandle.Layout.Tile = 'south';
end
end

function caseList = normalize_cases(cases)
if iscell(cases)
    caseList = cases;
elseif isstruct(cases)
    caseList = num2cell(cases);
else
    error('plot_prediction:InvalidCases', ...
        'predictionCases must be a cell array or struct array.');
end
end

function [numRows, numCols] = prediction_layout_size(numPanels)
if numPanels <= 2
    numRows = 1;
    numCols = max(1, numPanels);
else
    numCols = 2;
    numRows = ceil(numPanels / numCols);
end
end

function panels = build_prediction_panels(models, caseList)
panels = empty_panels();

for caseIndex = 1:length(caseList)
    caseData = caseList{caseIndex};
    seriesList = collect_case_series(models, caseData);

    for seriesIndex = 1:length(seriesList)
        series = seriesList(seriesIndex);
        panelKey = prediction_panel_key(series);
        panelIndex = find_panel(panels, panelKey);

        if panelIndex == 0
            panelIndex = length(panels) + 1;
            panels(panelIndex).key = panelKey;
            panels(panelIndex).component = series.component;
            panels(panelIndex).stress_label = stress_label(series.component);
            panels(panelIndex).seriesList = empty_series();
        end

        panels(panelIndex).seriesList(end + 1) = series;
    end
end
end

function panels = empty_panels()
panels = struct('key', {}, 'component', {}, 'stress_label', {}, ...
    'seriesList', {});
end

function key = prediction_panel_key(series)
key = upper(char(series.component));
end

function index = find_panel(panels, key)
index = 0;
for ii = 1:length(panels)
    if strcmpi(panels(ii).key, key)
        index = ii;
        return;
    end
end
end

function [experimentHandles, predictionHandles, experimentLabels, predictionLabels, legendKeys] = ...
    plot_panel(panel, axesHandle)
colors = scientific_colors();
markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};
experimentHandles = [];
predictionHandles = [];
experimentLabels = {};
predictionLabels = {};
legendKeys = {};

for seriesIndex = 1:length(panel.seriesList)
    color = colors(wrap_index(seriesIndex, size(colors, 1)), :);
    marker = markers{wrap_index(seriesIndex, length(markers))};
    series = panel.seriesList(seriesIndex);

    [experimentHandle, predictionHandle] = plot_prediction_series(axesHandle, ...
        series, color, marker);
    experimentHandles(end + 1) = experimentHandle;
    predictionHandles(end + 1) = predictionHandle;
    experimentLabels{end + 1} = sprintf('%s exp', series.legend_label);
    predictionLabels{end + 1} = sprintf('%s prediction', series.legend_label);
    legendKeys{end + 1} = series.legend_key;
end
end

function [sharedExperimentHandles, sharedPredictionHandles, ...
    sharedExperimentLabels, sharedPredictionLabels, sharedLegendKeys] = ...
    append_shared_legend_entries(sharedExperimentHandles, ...
    sharedPredictionHandles, sharedExperimentLabels, sharedPredictionLabels, ...
    sharedLegendKeys, experimentHandles, predictionHandles, ...
    experimentLabels, predictionLabels, legendKeys)
for ii = 1:length(legendKeys)
    if any(strcmp(sharedLegendKeys, legendKeys{ii}))
        continue;
    end

    sharedExperimentHandles(end + 1) = experimentHandles(ii);
    sharedPredictionHandles(end + 1) = predictionHandles(ii);
    sharedExperimentLabels{end + 1} = experimentLabels{ii};
    sharedPredictionLabels{end + 1} = predictionLabels{ii};
    sharedLegendKeys{end + 1} = legendKeys{ii};
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
seriesList(1).legend_label = compact_legend_label(caseData.name);
seriesList(1).legend_key = char(caseData.name);
seriesList(1).component = upper(char(caseData.stress_component));
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
    series.legend_label = compact_legend_label(caseData.name);
    series.legend_key = char(caseData.name);
    series.component = components{componentIndex};
    seriesList(end + 1) = series;
end
end

function seriesList = relative_cauchy_case_series(models, caseData)
seriesList = empty_series();
seriesList(1).lambda = caseData.lambda(:, 1);
seriesList(1).exp = caseData.sigma_exp(:);
seriesList(1).fit = fitted_relative_cauchy(models, caseData);
seriesList(1).label = sprintf('%s sigma11-sigma22', caseData.name);
seriesList(1).legend_label = compact_legend_label(caseData.name);
seriesList(1).legend_key = char(caseData.name);
seriesList(1).component = upper(char(caseData.relative_cauchy_component));
end

function seriesList = empty_series()
seriesList = struct('lambda', {}, 'exp', {}, 'fit', {}, 'label', {}, ...
    'legend_label', {}, 'legend_key', {}, 'component', {});
end

function label = compact_legend_label(name)
label = char(name);
token = regexp(label, 'lambda2=([0-9.]+)', 'tokens', 'once');
if ~isempty(token)
    label = sprintf('\\lambda_2 = %s', token{1});
end
end

function [experimentHandle, predictionHandle] = plot_prediction_series(axesHandle, ...
    series, color, marker)
[lambdaSorted, fitSorted] = sort_curve(series.lambda, series.fit);

experimentHandle = plot(axesHandle, series.lambda, series.exp, marker, ...
    'LineStyle', 'none', ...
    'MarkerSize', 5.5, ...
    'LineWidth', 1.1, ...
    'MarkerFaceColor', [1, 1, 1], ...
    'MarkerEdgeColor', color, ...
    'DisplayName', sprintf('%s exp', series.label));

predictionHandle = plot(axesHandle, lambdaSorted, fitSorted, '-', ...
    'Color', color, ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('%s prediction', series.label));
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
        error('plot_prediction:UnsupportedRelativeCauchyComponent', ...
            'Unsupported relative Cauchy component: %s.', name);
end
end

function label = stress_label(component)
switch upper(char(component))
    case 'P11'
        label = '$P_{11}$ (MPa)';
    case 'P22'
        label = '$P_{22}$ (MPa)';
    case 'P33'
        label = '$P_{33}$ (MPa)';
    case 'P12'
        label = '$P_{12}$ (MPa)';
    case 'P21'
        label = '$P_{21}$ (MPa)';
    case 'SIGMA11_MINUS_SIGMA22'
        label = '$\sigma_{11} - \sigma_{22}$ (MPa)';
    otherwise
        label = '$P$ (MPa)';
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
        error('plot_prediction:UnsupportedComponent', ...
            'Unsupported PK1 component: %s.', name);
end
end
