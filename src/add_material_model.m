function models = add_material_model(models, model, lowerBounds, upperBounds)
% Add one material model to a combined material model.
%
% Example:
%   models = [];
%   models = add_material_model(models, Yeoh([0.2, 0.01, 0.001]), ...
%       [0.0, -Inf, -Inf], [Inf, Inf, Inf]);
%   models = add_material_model(models, Neo_Hookean([0.1]), 0.0, Inf);

if nargin < 1 || isempty(models)
    terms = {};
else
    terms = models.terms;
end

if nargin < 2
    error('add_material_model:MissingModel', 'A material model must be provided.');
end

check_material_model(model);

if nargin < 3 || isempty(lowerBounds)
    if isfield(model, 'lower_bounds')
        lowerBounds = model.lower_bounds;
    else
        lowerBounds = -Inf(size(model.parameters));
    end
end
if nargin < 4 || isempty(upperBounds)
    if isfield(model, 'upper_bounds')
        upperBounds = model.upper_bounds;
    else
        upperBounds = Inf(size(model.parameters));
    end
end

lowerBounds = expand_bounds(lowerBounds, model.parameters, 'lowerBounds');
upperBounds = expand_bounds(upperBounds, model.parameters, 'upperBounds');

if any(lowerBounds > upperBounds)
    error('add_material_model:InvalidBounds', ...
        'Each lower bound must be smaller than or equal to the upper bound.');
end

model.lower_bounds = lowerBounds;
model.upper_bounds = upperBounds;

terms{end + 1} = model;
models = build_combined_model(terms);
end

function check_material_model(model)
requiredFields = {'name', 'parameters', 'parameter_names', 'energy', 'S', 'P', 'set_parameters'};
for ii = 1:length(requiredFields)
    fieldName = requiredFields{ii};
    if ~isfield(model, fieldName)
        error('add_material_model:InvalidModel', ...
            'Material model is missing field "%s".', fieldName);
    end
end

if length(model.parameter_names) ~= length(model.parameters)
    error('add_material_model:InvalidModel', ...
        'parameter_names must have the same length as parameters.');
end
end

function bounds = expand_bounds(bounds, parameters, name)
bounds = bounds(:).';
if isscalar(bounds)
    bounds = repmat(bounds, size(parameters));
elseif length(bounds) ~= length(parameters)
    error('add_material_model:InvalidBounds', ...
        '%s must be scalar or have length %d.', name, length(parameters));
end
end

function models = build_combined_model(terms)
labels = term_labels(terms);

models = struct();
models.name = combined_name(labels);
models.terms = terms;
models.parameters = collect_vector(terms, 'parameters');
models.lower_bounds = collect_vector(terms, 'lower_bounds');
models.upper_bounds = collect_vector(terms, 'upper_bounds');
models.parameter_names = collect_parameter_names(terms, labels);
models.parameter_map = collect_parameter_map(terms, labels);

models.energy = @(F) total_energy(terms, F);
models.S = @(F) total_S(terms, F);
models.P = @(F) incompressible_constraint(F * total_S(terms, F), F);
models.set_parameters = @(parameters) set_combined_parameters(terms, parameters);
models.objective = @(parameters, cases) objective(parameters, models, cases);
end

function labels = term_labels(terms)
baseNames = cell(1, length(terms));
for ii = 1:length(terms)
    baseNames{ii} = terms{ii}.name;
end

labels = baseNames;
for ii = 1:length(baseNames)
    repeated = sum(strcmp(baseNames, baseNames{ii})) > 1;
    if repeated
        occurrence = sum(strcmp(baseNames(1:ii), baseNames{ii}));
        labels{ii} = sprintf('%s(%d)', baseNames{ii}, occurrence);
    end
end
end

function name = combined_name(labels)
name = strjoin(labels, ' + ');
end

function values = collect_vector(terms, fieldName)
values = [];
for ii = 1:length(terms)
    values = [values, terms{ii}.(fieldName)];
end
end

function names = collect_parameter_names(terms, labels)
names = {};
for ii = 1:length(terms)
    term = terms{ii};
    for jj = 1:length(term.parameter_names)
        names{end + 1} = sprintf('%s.%s', labels{ii}, term.parameter_names{jj});
    end
end
end

function parameterMap = collect_parameter_map(terms, labels)
parameterMap = struct('index', {}, 'model', {}, 'parameter', {}, ...
    'initial_value', {}, 'lower_bound', {}, 'upper_bound', {});

index = 0;
for ii = 1:length(terms)
    term = terms{ii};
    for jj = 1:length(term.parameters)
        index = index + 1;
        parameterMap(index).index = index;
        parameterMap(index).model = labels{ii};
        parameterMap(index).parameter = term.parameter_names{jj};
        parameterMap(index).initial_value = term.parameters(jj);
        parameterMap(index).lower_bound = term.lower_bounds(jj);
        parameterMap(index).upper_bound = term.upper_bounds(jj);
    end
end
end

function W = total_energy(terms, F)
W = 0.0;
for ii = 1:length(terms)
    W = W + terms{ii}.energy(F);
end
end

function out = total_S(terms, F)
out = zeros(3, 3);
for ii = 1:length(terms)
    out = out + terms{ii}.S(F);
end
end

function models = set_combined_parameters(terms, parameters)
parameters = parameters(:).';
expectedLength = 0;
for ii = 1:length(terms)
    expectedLength = expectedLength + length(terms{ii}.parameters);
end

if length(parameters) ~= expectedLength
    error('add_material_model:InvalidParameterLength', ...
        'Expected %d parameters, got %d.', expectedLength, length(parameters));
end

newTerms = cell(size(terms));
startIndex = 1;
for ii = 1:length(terms)
    term = terms{ii};
    stopIndex = startIndex + length(term.parameters) - 1;
    termParameters = parameters(startIndex:stopIndex);

    newTerm = term.set_parameters(termParameters);
    newTerm.lower_bounds = term.lower_bounds;
    newTerm.upper_bounds = term.upper_bounds;

    newTerms{ii} = newTerm;
    startIndex = stopIndex + 1;
end

models = build_combined_model(newTerms);
end
