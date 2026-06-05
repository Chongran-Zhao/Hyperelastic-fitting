function textHandle = add_evaluation_text(axesHandle, metrics)
%   Place NMAD and R-square in the lower-right axes corner.
%
%   textHandle = ADD_EVALUATION_TEXT(axesHandle, metrics) adds a text label
%   to the specified axes showing the averaged fitting-quality indicators:
%
%       NMAD = normalized mean absolute difference
%       R^2  = coefficient of determination
%
%   The text is positioned using normalized axes coordinates, with
%   x = 0.98 and y = 0.04, so that it appears near the lower-right corner
%   independently of the current axis limits.
%
%   If metrics is empty or metrics.num_series is zero, no text is added and
%   an empty handle is returned.
%
%   The input metrics is expected to contain the fields
%
%       metrics.NMAD
%       metrics.R_square
%       metrics.num_series
%
%   as returned by evaluation_metrics.

textHandle = [];
if isempty(metrics) || metrics.num_series == 0
    return;
end

label = sprintf('NMAD = %.4g%%\nR^2 = %.4g', ...
    metrics.NMAD, metrics.R_square);

textHandle = text(axesHandle, 0.98, 0.04, label, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', ...
    'Interpreter', 'tex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'Color', [0.12, 0.12, 0.12]);
end