function [models, fit] = start_fit(models, cases, options)
% Fit one combined material model to all selected experiment cases.
%
% Example:
%   [models, fit] = start_fit(models, cases);

if nargin < 3 || isempty(options)
    options = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'FunctionTolerance', 1.0e-12, ...
        'StepTolerance', 1.0e-12, ...
        'MaxFunctionEvaluations', 2000);
end

iterationHistory = [];
resnormHistory = [];
iterationFigure = [];
iterationAxes = [];
iterationLine = [];

options = optimoptions(options, 'OutputFcn', @plot_iteration);

objectiveFunction = @(parameters) objective(parameters, models, cases);

[parameters, resnorm, residual, exitflag, output] = lsqnonlin( ...
    objectiveFunction, ...
    models.parameters, ...
    models.lower_bounds, ...
    models.upper_bounds, ...
    options);

models = models.set_parameters(parameters);

fit = struct();
fit.parameters = parameters;
fit.resnorm = resnorm;
fit.residual = residual;
fit.exitflag = exitflag;
fit.output = output;
fit.options = options;
fit.iterations = iterationHistory;
fit.resnorm_history = resnormHistory;

print_fitted_parameters(models);

    function stop = plot_iteration(parameters, optimValues, state)
        %#ok<INUSD>
        stop = false;

        switch state
            case 'init'
                initialize_iteration_plot();
            case 'iter'
                append_iteration_point(optimValues);
            case 'done'
                append_iteration_point(optimValues);
        end
    end

    function initialize_iteration_plot()
        iterationHistory = [];
        resnormHistory = [];

        iterationFigure = figure('Name', 'Fitting iteration', 'Color', 'w');
        iterationAxes = axes(iterationFigure);
        hold(iterationAxes, 'on');
        iterationLine = plot(iterationAxes, NaN, NaN, '-o', ...
            'Color', [0.00, 0.45, 0.74], ...
            'MarkerFaceColor', [1, 1, 1], ...
            'LineWidth', 1.8, ...
            'MarkerSize', 5);

        xlabel(iterationAxes, 'Iteration', ...
            'FontName', 'Times New Roman');
        ylabel(iterationAxes, 'Residual norm', ...
            'FontName', 'Times New Roman');
        set(iterationAxes, ...
            'Box', 'off', ...
            'LineWidth', 1.1, ...
            'FontName', 'Times New Roman', ...
            'FontSize', 12, ...
            'TickDir', 'out', ...
            'YScale', 'log');
        grid(iterationAxes, 'on');
    end

    function append_iteration_point(optimValues)
        if isempty(iterationAxes) || ~isvalid(iterationAxes)
            initialize_iteration_plot();
        end

        if ~isfield(optimValues, 'iteration') || ...
                ~isfield(optimValues, 'resnorm')
            return;
        end

        iterationHistory(end + 1) = optimValues.iteration;
        resnormHistory(end + 1) = max(optimValues.resnorm, realmin);

        if length(iterationHistory) >= 2 && ...
                iterationHistory(end) == iterationHistory(end - 1)
            iterationHistory(end - 1) = [];
            resnormHistory(end - 1) = [];
        end

        set(iterationLine, ...
            'XData', iterationHistory, ...
            'YData', resnormHistory);
        drawnow limitrate;
    end
end

function print_fitted_parameters(models)
fprintf('\nOptimized parameters:\n');

parameterMap = models.parameter_map;
if isempty(parameterMap)
    return;
end

currentModel = '';
branchIndex = 0;

for ii = 1:length(parameterMap)
    item = parameterMap(ii);

    if ~strcmp(currentModel, item.model)
        branchIndex = branchIndex + 1;
        currentModel = item.model;
        fprintf('branch %d: %s\n', branchIndex, currentModel);
    end

    fprintf('  %s = %.12g\n', item.parameter, item.initial_value);
end
end
