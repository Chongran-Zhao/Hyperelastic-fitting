function cases = add_exp_data_sets(cases, author, mode, varargin)
%add_exp_data_sets Append a registered experiment case.
%
% Examples:
%   cases = add_exp_data_sets(cases, 'Treloar', 'UT');
%   cases = add_exp_data_sets(cases, 'Kawamura', 'BT', 1.7);
%
% Call add_exp_data_sets() without inputs to print registered calls.

if nargin == 0
    cases = registered_data_calls();
    if nargout == 0
        print_registered_data_calls(cases);
    end
    return;
end

if isempty(cases)
    cases = {};
end

if nargin < 3
    error('add_exp_data_sets:MissingInput', ...
        'Use add_exp_data_sets(cases, author, mode) or add_exp_data_sets(cases, author, ''BT'', fixedStretch).');
end

spec = data_case_spec(author, mode, varargin{:});
cases{end + 1} = load_experiment_case(spec);
end

function calls = registered_data_calls()
calls = { ...
    'add_exp_data_sets(cases, ''Treloar'', ''UT'')', ...
    'add_exp_data_sets(cases, ''Treloar'', ''ET'')', ...
    'add_exp_data_sets(cases, ''Treloar'', ''PS'')', ...
    'add_exp_data_sets(cases, ''Kawabata'', ''UT'')', ...
    'add_exp_data_sets(cases, ''Kawabata'', ''ET'')', ...
    'add_exp_data_sets(cases, ''Kawabata'', ''PS'')', ...
    'add_exp_data_sets(cases, ''Meunier'', ''UT'')', ...
    'add_exp_data_sets(cases, ''Meunier'', ''ET'')', ...
    'add_exp_data_sets(cases, ''Meunier'', ''PS'')', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''UT'')', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''UE'')', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''ET'')', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''BT'', 1.7)', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''BT'', 1.5)', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''BT'', 1.3)', ...
    'add_exp_data_sets(cases, ''Kawamura'', ''BT'', 1.1)', ...
    'add_exp_data_sets(cases, ''Jones'', ''UT'')', ...
    'add_exp_data_sets(cases, ''Jones'', ''BT'', 1.0)', ...
    'add_exp_data_sets(cases, ''Jones'', ''BT'', 1.502)', ...
    'add_exp_data_sets(cases, ''Jones'', ''BT'', 1.984)', ...
    'add_exp_data_sets(cases, ''Jones'', ''BT'', 2.295)', ...
    'add_exp_data_sets(cases, ''Jones'', ''BT'', 2.623)', ...
    'add_exp_data_sets(cases, ''James'', ''UT'')', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 1.3)', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 1.5)', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 1.7)', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 2.0)', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 2.5)', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 3.0)', ...
    'add_exp_data_sets(cases, ''James'', ''BT'', 3.5)', ...
    'add_exp_data_sets(cases, ''Katashima'', ''ET'')', ...
    'add_exp_data_sets(cases, ''Katashima'', ''BT'')', ...
    'add_exp_data_sets(cases, ''Katashima'', ''PS'')' ...
};
end

function print_registered_data_calls(calls)
for ii = 1:length(calls)
    fprintf('%s\n', calls{ii});
end
end

function spec = data_case_spec(author, mode, varargin)
authorId = normalize_name(author);
modeId = upper(char(mode));
rootDir = fileparts(fileparts(mfilename('fullpath')));

switch authorId
    case 'treloar'
        spec = table_family_spec('Treloar', modeId, rootDir, 'Treloar_1944');
    case 'kawabata'
        spec = table_family_spec('Kawabata', modeId, rootDir, 'Kawabata_1981');
    case 'meunier'
        spec = table_family_spec('Meunier', modeId, rootDir, 'Meunier_2008');
    case 'kawamura'
        spec = kawamura_spec(modeId, rootDir, varargin{:});
    case 'jones'
        spec = jones_spec(modeId, rootDir, varargin{:});
    case 'james'
        spec = james_spec(modeId, rootDir, varargin{:});
    case 'katashima'
        spec = katashima_spec(modeId, rootDir, varargin{:});
    otherwise
        error('add_exp_data_sets:UnknownAuthor', ...
            'Unknown data author: %s. Call add_exp_data_sets() to see registered cases.', ...
            char(author));
end
end

function id = normalize_name(value)
id = lower(char(value));
id = strrep(id, '-', '_');
id = strrep(id, ' ', '_');
end

function spec = table_family_spec(name, mode, rootDir, folderName)
switch mode
    case {'UT', 'ET', 'PS'}
        spec = table_spec(sprintf('%s %s', name, mode), mode, rootDir, ...
            sprintf('%s/%s/stress_stretch.txt', folderName, mode));
    otherwise
        error('add_exp_data_sets:UnsupportedMode', ...
            '%s data supports UT, ET, and PS; got %s.', name, mode);
end
end

function spec = kawamura_spec(mode, rootDir, varargin)
switch mode
    case {'UT', 'UE', 'ET'}
        spec = vector_spec(sprintf('Kawamura %s', mode), mode, rootDir, ...
            sprintf('Kawamura_2001/%s/stretch.txt', mode), ...
            sprintf('Kawamura_2001/%s/stress.txt', mode), 'P11');
    case 'BT'
        if isempty(varargin)
            error('add_exp_data_sets:MissingFixedStretch', ...
                'Use add_exp_data_sets(cases, ''Kawamura'', ''BT'', fixedStretch).');
        end
        spec = kawamura_bt_spec(rootDir, varargin{1});
    otherwise
        error('add_exp_data_sets:UnsupportedMode', ...
            'Kawamura data supports UT, UE, ET, and BT; got %s.', mode);
end
end

function spec = jones_spec(mode, rootDir, varargin)
switch mode
    case 'UT'
        spec = vector_spec('Jones UT', 'UT', rootDir, ...
            'Jones_1975/UT/stretch.txt', ...
            'Jones_1975/UT/stress.txt', 'P11');
    case 'BT'
        if isempty(varargin)
            error('add_exp_data_sets:MissingFixedStretch', ...
                'Use add_exp_data_sets(cases, ''Jones'', ''BT'', fixedStretch).');
        end
        spec = jones_bt_spec(rootDir, varargin{1});
    otherwise
        error('add_exp_data_sets:UnsupportedMode', ...
            'Jones data supports UT and BT; got %s.', mode);
end
end

function spec = james_spec(mode, rootDir, varargin)
switch mode
    case 'UT'
        spec = vector_spec('James UT', 'UT', rootDir, ...
            'James_1975/UT/stretch.txt', ...
            'James_1975/UT/stress.txt', 'P11');
    case 'BT'
        if isempty(varargin)
            error('add_exp_data_sets:MissingFixedStretch', ...
                'Use add_exp_data_sets(cases, ''James'', ''BT'', fixedStretch).');
        end
        spec = james_bt_spec(rootDir, varargin{1});
    otherwise
        error('add_exp_data_sets:UnsupportedMode', ...
            'James data supports UT and BT; got %s.', mode);
end
end

function spec = katashima_spec(mode, rootDir, varargin)
switch mode
    case 'ET'
        spec = vector_spec('Katashima ET', 'ET', rootDir, ...
            'Katashima_2012/ET/stretch.txt', ...
            'Katashima_2012/ET/stress.txt', 'P11');
    case 'BT'
        spec = katashima_bt_spec(rootDir);
    case 'PS'
        spec = katashima_ps_spec(rootDir);
    otherwise
        error('add_exp_data_sets:UnsupportedMode', ...
            'Katashima data supports ET, BT, and PS; got %s.', mode);
end
end

function spec = table_spec(name, mode, rootDir, relPath)
spec = struct();
spec.name = name;
spec.source = relPath;
spec.mode = mode;
spec.file = data_path(rootDir, relPath);
spec.columns.P11 = 1;
spec.columns.lambda1 = 2;
spec.stressComponent = 'P11';
end

function spec = vector_spec(name, mode, rootDir, stretchRelPath, stressRelPath, component)
component = upper(component);
spec = struct();
spec.name = name;
spec.source = stressRelPath;
spec.mode = mode;
spec.lambda1File = data_path(rootDir, stretchRelPath);
spec.(sprintf('%sFile', component)) = data_path(rootDir, stressRelPath);
spec.stressComponent = component;
end

function spec = kawamura_bt_spec(rootDir, fixedStretch)
[lambdaLabel, lambda2] = fixed_stretch_label(fixedStretch);
baseRelPath = sprintf('Kawamura_2001/BT/lambda_%s', lambdaLabel);

obs(1).component = 'P11';
obs(1).lambda1File = data_path(rootDir, sprintf('%s/stretch_1.txt', baseRelPath));
obs(1).stressFile = data_path(rootDir, sprintf('%s/stress_1.txt', baseRelPath));
obs(1).lambda2 = lambda2;

obs(2).component = 'P22';
obs(2).lambda1File = data_path(rootDir, sprintf('%s/stretch_1_2.txt', baseRelPath));
obs(2).stressFile = data_path(rootDir, sprintf('%s/stress_2.txt', baseRelPath));
obs(2).lambda2 = lambda2;

spec = struct();
spec.name = sprintf('Kawamura BT lambda2=%g', lambda2);
spec.source = baseRelPath;
spec.mode = 'BT';
spec.lambda2 = lambda2;
spec.observations = obs;
end

function spec = jones_bt_spec(rootDir, fixedStretch)
knownValues = [1.0, 1.502, 1.984, 2.295, 2.623];
knownLabels = {'1', '1d502', '1d984', '2d295', '2d623'};
[lambdaLabel, lambda2] = fixed_stretch_label( ...
    fixedStretch, knownValues, knownLabels, 'Jones BT');

baseRelPath = sprintf('Jones_1975/Biaxial_tension/lambda_%s', lambdaLabel);

spec = relative_cauchy_spec( ...
    sprintf('Jones BT lambda2=%g', lambda2), ...
    'BT', rootDir, ...
    sprintf('%s/stretch.txt', baseRelPath), ...
    sprintf('%s/relative_stress.txt', baseRelPath), ...
    lambda2, ...
    'SIGMA11_MINUS_SIGMA22');
end

function spec = james_bt_spec(rootDir, fixedStretch)
knownValues = [1.3, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5];
knownLabels = {'1d3', '1d5', '1d7', '2d0', '2d5', '3d0', '3d5'};
[lambdaLabel, lambda2] = fixed_stretch_label( ...
    fixedStretch, knownValues, knownLabels, 'James BT');

baseRelPath = sprintf('James_1975/BT/lambda_%s', lambdaLabel);

obs(1).component = 'P11';
obs(1).lambda1File = data_path(rootDir, sprintf('%s/stretch_1.txt', baseRelPath));
obs(1).stressFile = data_path(rootDir, sprintf('%s/stress_1.txt', baseRelPath));
obs(1).lambda2 = lambda2;

obs(2).component = 'P22';
obs(2).lambda1File = data_path(rootDir, sprintf('%s/stretch_1_2.txt', baseRelPath));
obs(2).stressFile = data_path(rootDir, sprintf('%s/stress_2.txt', baseRelPath));
obs(2).lambda2 = lambda2;

spec = struct();
spec.name = sprintf('James BT lambda2=%g', lambda2);
spec.source = baseRelPath;
spec.mode = 'BT';
spec.lambda2 = lambda2;
spec.observations = obs;
end

function spec = katashima_bt_spec(rootDir)
baseRelPath = 'Katashima_2012/BT';

obs(1).component = 'P11';
obs(1).lambda1File = data_path(rootDir, sprintf('%s/stretch_1.txt', baseRelPath));
obs(1).stressFile = data_path(rootDir, sprintf('%s/stress_1.txt', baseRelPath));
obs(1).lambda2 = @(lambda1) 0.5 + 0.5 .* lambda1;

obs(2).component = 'P22';
obs(2).lambda1File = data_path(rootDir, sprintf('%s/stretch_1_2.txt', baseRelPath));
obs(2).stressFile = data_path(rootDir, sprintf('%s/stress_2.txt', baseRelPath));
obs(2).lambda2 = @(lambda1) 0.5 + 0.5 .* lambda1;

spec = struct();
spec.name = 'Katashima BT';
spec.source = baseRelPath;
spec.mode = 'BT';
spec.observations = obs;
end

function spec = katashima_ps_spec(rootDir)
baseRelPath = 'Katashima_2012/PS';

obs(1).component = 'P11';
obs(1).lambda1File = data_path(rootDir, sprintf('%s/stretch_1.txt', baseRelPath));
obs(1).stressFile = data_path(rootDir, sprintf('%s/stress_1.txt', baseRelPath));

obs(2).component = 'P22';
obs(2).lambda1File = data_path(rootDir, sprintf('%s/stretch_1_2.txt', baseRelPath));
obs(2).stressFile = data_path(rootDir, sprintf('%s/stress_2.txt', baseRelPath));

spec = struct();
spec.name = 'Katashima PS';
spec.source = baseRelPath;
spec.mode = 'PS';
spec.observations = obs;
end

function spec = relative_cauchy_spec(name, mode, rootDir, stretchRelPath, ...
    stressRelPath, lambda2, component)
spec = struct();
spec.name = name;
spec.source = stressRelPath;
spec.mode = mode;
spec.lambda1File = data_path(rootDir, stretchRelPath);
spec.relativeCauchyFile = data_path(rootDir, stressRelPath);
spec.relativeCauchyComponent = component;
spec.lambda2 = lambda2;
end

function [label, value] = fixed_stretch_label(fixedStretch, knownValues, knownLabels, dataName)
value = fixedStretch;
if ischar(value)
    value = str2double(value);
elseif (exist('isstring', 'builtin') || exist('isstring', 'file')) && isstring(value)
    value = str2double(char(value));
end

if nargin < 2
    knownValues = [1.7, 1.5, 1.3, 1.1];
    knownLabels = {'1d7', '1d5', '1d3', '1d1'};
    dataName = 'Kawamura BT';
end

[difference, index] = min(abs(knownValues - value));
if difference > 1.0e-12
    error('add_exp_data_sets:UnsupportedFixedStretch', ...
        '%s does not support fixed stretch %g.', dataName, value);
end

label = knownLabels{index};
value = knownValues(index);
end

function caseData = load_experiment_case(spec)
mode = char(require_field(spec, 'mode'));
stressScale = get_optional_field(spec, 'stressScale', 1.0);

if isfield(spec, 'relativeCauchyFile')
    caseData = load_relative_cauchy_case(spec, stressScale);
    return;
end

if isfield(spec, 'observations')
    caseData = load_observation_case(spec, stressScale);
    return;
end

tableData = [];
if isfield(spec, 'file')
    tableData = read_numeric_file(spec.file);
end

lambda1 = read_signal(spec, tableData, 'lambda1', ...
    {'lambda1File', 'stretchFile'}, []);
if isempty(lambda1)
    error('add_exp_data_sets:MissingLambda1', 'spec must provide lambda1 data.');
end
lambda1 = lambda1(:);

lambda2 = read_signal(spec, tableData, 'lambda2', {'lambda2File'}, []);
if isempty(lambda2) && isfield(spec, 'lambda2')
    lambda2 = spec.lambda2;
end

F_list = build_F_list(mode, lambda1, lambda2);
n = length(lambda1);

P11 = read_signal(spec, tableData, 'P11', {'P11File', 'stressFile'}, []);
P22 = read_signal(spec, tableData, 'P22', {'P22File'}, []);
P33 = read_signal(spec, tableData, 'P33', {'P33File'}, []);
stressComponent = char(get_optional_field(spec, 'stressComponent', 'P11'));

if ~isempty(P11)
    P11 = stressScale .* P11;
end
if ~isempty(P22)
    P22 = stressScale .* P22;
end
if ~isempty(P33)
    P33 = stressScale .* P33;
end

caseData = struct();
caseData.name = get_optional_field(spec, 'name', '');
caseData.source = get_optional_field(spec, 'source', '');
caseData.mode = mode;
caseData.lambda = get_principal_stretches(F_list);
caseData.F_list = F_list;
caseData.num_points = n;

if strcmpi(mode, 'BT')
    [P_list, P_mask] = build_PK1_list(n, ...
        'P11', P11, ...
        'P22', P22, ...
        'P33', P33);
    caseData.stress_type = 'tensor';
    caseData.P_list = P_list;
    caseData.P_mask = P_mask;
else
    caseData.stress_type = 'scalar';
    caseData.stress_component = stressComponent;
    caseData.P_exp = get_scalar_stress(stressComponent, P11, P22, P33);
end
end

function caseData = load_observation_case(spec, stressScale)
mode = char(spec.mode);
F_list = zeros(3, 3, 0);
P_list = zeros(3, 3, 0);
P_mask = false(3, 3, 0);

for ii = 1:length(spec.observations)
    observation = spec.observations(ii);
    lambda1 = read_numeric_file(observation.lambda1File);
    lambda1 = lambda1(:);

    if isfield(observation, 'lambda2')
        lambda2 = evaluate_lambda2(observation.lambda2, lambda1);
    elseif isfield(spec, 'lambda2')
        lambda2 = evaluate_lambda2(spec.lambda2, lambda1);
    else
        lambda2 = [];
    end

    stress = read_numeric_file(observation.stressFile);
    stress = stressScale .* stress(:);

    F_obs = build_F_list(mode, lambda1, lambda2);
    [P_obs, mask_obs] = build_PK1_list(length(lambda1), ...
        observation.component, stress);

    F_list = cat(3, F_list, F_obs);
    P_list = cat(3, P_list, P_obs);
    P_mask = cat(3, P_mask, mask_obs);
end

caseData = struct();
caseData.name = get_optional_field(spec, 'name', '');
caseData.source = get_optional_field(spec, 'source', '');
caseData.mode = mode;
caseData.lambda = get_principal_stretches(F_list);
caseData.F_list = F_list;
caseData.stress_type = 'tensor';
caseData.P_list = P_list;
caseData.P_mask = P_mask;
caseData.num_points = size(F_list, 3);
end

function caseData = load_relative_cauchy_case(spec, stressScale)
mode = char(spec.mode);
lambda1 = read_numeric_file(spec.lambda1File);
lambda1 = lambda1(:);
lambda2 = evaluate_lambda2(spec.lambda2, lambda1);
F_list = build_F_list(mode, lambda1, lambda2);

sigmaDifference = read_numeric_file(spec.relativeCauchyFile);
sigmaDifference = stressScale .* sigmaDifference(:);

caseData = struct();
caseData.name = get_optional_field(spec, 'name', '');
caseData.source = get_optional_field(spec, 'source', '');
caseData.mode = mode;
caseData.lambda = get_principal_stretches(F_list);
caseData.F_list = F_list;
caseData.stress_type = 'relative_cauchy';
caseData.relative_cauchy_component = spec.relativeCauchyComponent;
caseData.sigma_exp = sigmaDifference;
caseData.num_points = size(F_list, 3);
end

function lambda2 = evaluate_lambda2(lambda2Spec, lambda1)
if isa(lambda2Spec, 'function_handle')
    lambda2 = lambda2Spec(lambda1);
else
    lambda2 = lambda2Spec;
end
end

function F_list = build_F_list(mode, lambda1, lambda2)
mode = char(mode);
lambda1 = lambda1(:);

switch upper(mode)
    case {'UT', 'UE'}
        lambda2 = lambda1.^(-0.5);
        lambda3 = lambda1.^(-0.5);
    case 'ET'
        lambda2 = lambda1;
        lambda3 = lambda1.^(-2.0);
    case 'PS'
        lambda2 = ones(size(lambda1));
        lambda3 = lambda1.^(-1.0);
    case 'BT'
        if nargin < 3
            error('add_exp_data_sets:MissingLambda2', 'BT mode requires lambda2.');
        end
        lambda2 = expand_to_length(lambda2, length(lambda1), 'lambda2');
        lambda3 = 1.0 ./ (lambda1 .* lambda2);
    otherwise
        error('add_exp_data_sets:UnsupportedMode', ...
            'Unsupported loading mode: %s.', mode);
end

n = length(lambda1);
F_list = zeros(3, 3, n);
F_list(1, 1, :) = lambda1;
F_list(2, 2, :) = lambda2;
F_list(3, 3, :) = lambda3;
end

function [P_list, P_mask] = build_PK1_list(n, varargin)
P_list = zeros(3, 3, n);
P_mask = false(3, 3, n);

if mod(length(varargin), 2) ~= 0
    error('add_exp_data_sets:InvalidStressArguments', ...
        'Stress components must be provided as name/value pairs.');
end

for ii = 1:2:length(varargin)
    name = varargin{ii};
    value = varargin{ii + 1};

    if isempty(value)
        continue;
    end

    values = expand_to_length(value, n, name);
    [row, col] = component_indices(name);
    P_list(row, col, :) = values;
    P_mask(row, col, :) = true;
end
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
        error('add_exp_data_sets:UnsupportedComponent', ...
            'Unsupported PK1 component: %s.', name);
end
end

function value = read_signal(spec, tableData, columnName, fileFields, defaultValue)
value = defaultValue;

if isfield(spec, columnName)
    value = spec.(columnName);
    return;
end

if ~isempty(tableData) && isfield(spec, 'columns') && ...
        isfield(spec.columns, columnName)
    columnIndex = spec.columns.(columnName);
    value = tableData(:, columnIndex);
    return;
end

for ii = 1:length(fileFields)
    fieldName = fileFields{ii};
    if isfield(spec, fieldName)
        value = read_numeric_file(spec.(fieldName));
        value = value(:);
        return;
    end
end
end

function data = read_numeric_file(filePath)
filePath = char(filePath);
data = importdata(filePath);
if isstruct(data)
    data = data.data;
end
if isempty(data)
    error('add_exp_data_sets:EmptyFile', ...
        'No numeric data found in %s.', filePath);
end
end

function stress = get_scalar_stress(component, P11, P22, P33)
switch upper(component)
    case 'P11'
        stress = P11;
    case 'P22'
        stress = P22;
    case 'P33'
        stress = P33;
    otherwise
        error('add_exp_data_sets:UnsupportedScalarComponent', ...
            'Unsupported scalar PK1 component: %s.', component);
end

if isempty(stress)
    error('add_exp_data_sets:MissingStress', ...
        'No scalar stress data found for component %s.', component);
end

stress = stress(:);
end

function out = expand_to_length(value, n, name)
out = value(:);
if isscalar(out)
    out = repmat(out, n, 1);
elseif length(out) ~= n
    error('add_exp_data_sets:LengthMismatch', ...
        '%s must have length 1 or %d, got %d.', name, n, length(out));
end
end

function value = require_field(spec, name)
if ~isfield(spec, name)
    error('add_exp_data_sets:MissingField', 'spec.%s is required.', name);
end
value = spec.(name);
end

function value = get_optional_field(spec, name, defaultValue)
if isfield(spec, name)
    value = spec.(name);
else
    value = defaultValue;
end
end

function lambda = get_principal_stretches(F_list)
n = size(F_list, 3);
lambda = zeros(n, 3);
lambda(:, 1) = reshape(F_list(1, 1, :), n, 1);
lambda(:, 2) = reshape(F_list(2, 2, :), n, 1);
lambda(:, 3) = reshape(F_list(3, 3, :), n, 1);
end

function out = data_path(rootDir, relPath)
out = fullfile(rootDir, 'data', relPath);
end
