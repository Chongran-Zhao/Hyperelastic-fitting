function model = Micro_GenStrain(varargin)
% General chain-network model based on generalized strains.
%
% The chain statistics, chain-level strain E_hat(lambda_n), and
% macroscopic generalized strain E_bar are independently selectable. For
% each quadrature direction n,
%
%   E_hat(lambda_n) = E_bar:(n*n')
%
% and lambda_n is obtained from the inverse of E_hat.
%
% The unified stress formula is
%
%   S = J^(-2/3) DEV(S_bar),
%   S_bar = dPsi/dE_bar : Q_bar,
%
% where Q_bar = 2*dE_bar/dC_bar. The sphere average is evaluated with the
% Lebedev quadrature.
%
% Supported chain statistics:
%   Gaussian, NonGaussian.
%
% Supported chain and macroscopic generalized strains:
%   SH, Hencky, Biot, BI, CR, CZ, DN.
%
% Examples:
%   Micro_GenStrain('Gaussian', 'Hencky', 'Hencky')
%   Micro_GenStrain('NonGaussian', 'BI', 'Biot')
%   Micro_GenStrain([1.0, 50.0, 1.0], 'NonGaussian', 'BI', 'Biot')
%
% Explicit parameter vectors:
%   Gaussian:    parameters = [mu, chain parameters..., strain parameters...]
%   NonGaussian: parameters = [mu, N, chain parameters..., strain parameters...]

[parameters, statisticsFamily, chainFamily, strainFamily, ...
    lowerBounds, upperBounds] = parse_inputs(varargin{:});

statisticsParameterCount = statistics_parameter_count(statisticsFamily);
chainParameterCount = strain_parameter_count(chainFamily);
strainParameterCount = strain_parameter_count(strainFamily);
expectedParameterCount = statisticsParameterCount + ...
    chainParameterCount + strainParameterCount;

if length(parameters) ~= expectedParameterCount
    error('Micro_GenStrain:InvalidParameters', ...
        'Expected %d parameters, got %d.', ...
        expectedParameterCount, length(parameters));
end

[statisticsParameters, chainParameters, strainParameters] = ...
    split_parameters(parameters, statisticsParameterCount, ...
    chainParameterCount, strainParameterCount);

statisticsInfo = chain_statistics(statisticsFamily, statisticsParameters);
chainInfo = generalized_strain(chainFamily, ones(3, 1), chainParameters);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('MicroGenStrain-%s-Ehat%s-E%s', ...
    statisticsInfo.family, chainInfo.family, strainInfo.family);
model.parameters = parameters;
model.lower_bounds = lowerBounds;
model.upper_bounds = upperBounds;
model.statistics_family = statisticsInfo.family;
model.statistics_family_label = statisticsInfo.family_label;
model.chain_family = chainInfo.family;
model.chain_family_label = chainInfo.family_label;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [statisticsInfo.parameter_names, ...
    prefixed_parameter_names('chain', chainInfo.parameter_names), ...
    prefixed_parameter_names('strain', strainInfo.parameter_names)];
model.energy = @(F) Energy(parameters, statisticsFamily, ...
    chainFamily, strainFamily, F);
model.S = @(F) S(parameters, statisticsFamily, ...
    chainFamily, strainFamily, F);
model.P = @(F) P(parameters, statisticsFamily, ...
    chainFamily, strainFamily, F);
model.set_parameters = @(parameters) Micro_GenStrain( ...
    parameters, statisticsFamily, chainFamily, strainFamily);
end

function [parameters, statisticsFamily, chainFamily, strainFamily, ...
    lowerBounds, upperBounds] = parse_inputs(varargin)
inputCount = nargin;

if inputCount == 0
    error('Micro_GenStrain:InvalidCall', ...
        ['Use Micro_GenStrain(statisticsFamily, chainFamily, ', ...
        'strainFamily) for default parameters.']);
end

firstInput = varargin{1};

if is_text(firstInput)
    if inputCount ~= 3
        error('Micro_GenStrain:InvalidCall', ...
            ['Use Micro_GenStrain(statisticsFamily, chainFamily, ', ...
            'strainFamily) for default parameters.']);
    end

    statisticsFamily = firstInput;
    chainFamily = varargin{2};
    strainFamily = varargin{3};
    [parameters, lowerBounds, upperBounds] = default_parameters( ...
        statisticsFamily, chainFamily, strainFamily);
    return;
end

if inputCount ~= 4
    error('Micro_GenStrain:InvalidCall', ...
        ['Use Micro_GenStrain(parameters, statisticsFamily, chainFamily, ', ...
        'strainFamily) for explicit parameters.']);
end

parameters = firstInput(:).';
statisticsFamily = varargin{2};
chainFamily = varargin{3};
strainFamily = varargin{4};
[~, lowerBounds, upperBounds] = default_parameters( ...
    statisticsFamily, chainFamily, strainFamily);
end

function out = is_text(value)
out = ischar(value) || ...
    ((exist('isstring', 'builtin') || exist('isstring', 'file')) && isstring(value));
end

function [parameters, lowerBounds, upperBounds] = default_parameters( ...
    statisticsFamily, chainFamily, strainFamily)
[statisticsParameters, statisticsLowerBounds, statisticsUpperBounds] = ...
    default_statistics_parameters(statisticsFamily);
[chainParameters, chainLowerBounds, chainUpperBounds] = ...
    default_family_parameters(chainFamily);
[strainParameters, strainLowerBounds, strainUpperBounds] = ...
    default_family_parameters(strainFamily);

parameters = [statisticsParameters, chainParameters, strainParameters];
lowerBounds = [statisticsLowerBounds, chainLowerBounds, strainLowerBounds];
upperBounds = [statisticsUpperBounds, chainUpperBounds, strainUpperBounds];
end

function [parameters, lowerBounds, upperBounds] = ...
    default_statistics_parameters(statisticsFamily)
statisticsId = normalize_statistics_family(statisticsFamily);

switch statisticsId
    case 'gaussian'
        parameters = 1.0;
        lowerBounds = 0.0;
        upperBounds = Inf;

    case 'nongaussian'
        parameters = [1.0, 50.0];
        lowerBounds = [0.0, 1.0 + 1.0e-8];
        upperBounds = [Inf, Inf];

    otherwise
        error('Micro_GenStrain:UnsupportedStatistics', ...
            'Unsupported chain statistics: %s.', char(statisticsFamily));
end
end

function [parameters, lowerBounds, upperBounds] = ...
    default_family_parameters(family)
familyId = normalize_family(family);

switch familyId
    case {'hencky', 'biot'}
        parameters = [];
        lowerBounds = [];
        upperBounds = [];

    case {'sh', 'bi'}
        parameters = 1.0;
        lowerBounds = -Inf;
        upperBounds = Inf;

    case 'cr'
        parameters = [1.0, 1.0];
        lowerBounds = [1.0e-8, 1.0e-8];
        upperBounds = [Inf, Inf];

    case 'cz'
        parameters = 0.0;
        lowerBounds = -2.0;
        upperBounds = 2.0;

    case 'dn'
        parameters = [1.0, 1.0];
        lowerBounds = [1.0e-8, 1.0e-8];
        upperBounds = [Inf, Inf];

    otherwise
        error('Micro_GenStrain:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
end
end

function W = Energy(parameters, statisticsFamily, chainFamily, strainFamily, F)
[statisticsParameters, chainParameters, strainParameters] = ...
    split_parameters_by_family(parameters, statisticsFamily, ...
    chainFamily, strainFamily);

statistics = chain_statistics(statisticsFamily, statisticsParameters);
kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
chain = generalized_strain(chainFamily, ones(3, 1), chainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

W0 = chain_energy_density(1.0, statistics);
energyDensity = @(theta, phi) ...
    chain_energy_density(chain_stretch(E_bar, chain, theta, phi), ...
    statistics) - W0;
W = lebedev_quadrature(energyDensity);
end

function out = S(parameters, statisticsFamily, chainFamily, strainFamily, F)
[statisticsParameters, chainParameters, strainParameters] = ...
    split_parameters_by_family(parameters, statisticsFamily, ...
    chainFamily, strainFamily);

statistics = chain_statistics(statisticsFamily, statisticsParameters);
C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);
chain = generalized_strain(chainFamily, ones(3, 1), chainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

dPsi_dE_bar = lebedev_quadrature( ...
    @(theta, phi) chain_stress_tensor( ...
    E_bar, chain, chainFamily, chainParameters, statistics, theta, phi));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function tensor = chain_stress_tensor(E_bar, chain, chainFamily, ...
    chainParameters, statistics, theta, phi)
lambda_n = chain_stretch(E_bar, chain, theta, phi);
if ~is_admissible_stretch(lambda_n, statistics)
    tensor = NaN .* zeros(3, 3);
    return;
end

chainAtStretch = generalized_strain(chainFamily, lambda_n, chainParameters);
dW_dlambda = chain_energy_derivative(lambda_n, statistics);
g_lambda = 1.0 ./ chainAtStretch.derivatives(1);
tensor = dW_dlambda .* g_lambda .* direction_tensor(theta, phi);
end

function lambda_n = chain_stretch(E_bar, chain, theta, phi)
xi_n = contract(E_bar, direction_tensor(theta, phi));
lambda_n = chain.inverse(xi_n);
lambda_n = lambda_n(1);

if ~isreal(lambda_n) || lambda_n <= 0.0
    lambda_n = NaN;
end
end

function W = chain_energy_density(lambda_n, statistics)
if ~is_admissible_stretch(lambda_n, statistics)
    W = NaN;
    return;
end

switch statistics.family_id
    case 'gaussian'
        mu = statistics.parameters(1);
        W = 1.5 .* mu .* lambda_n .^ 2.0;

    case 'nongaussian'
        mu = statistics.parameters(1);
        N = statistics.parameters(2);
        beta = inv_langevin(lambda_n ./ sqrt(N));
        W = mu .* N .* ...
            (lambda_n .* beta ./ sqrt(N) + log(beta ./ sinh(beta)));

    otherwise
        error('Micro_GenStrain:UnsupportedStatistics', ...
            'Unsupported chain statistics: %s.', statistics.family);
end
end

function dW = chain_energy_derivative(lambda_n, statistics)
if ~is_admissible_stretch(lambda_n, statistics)
    dW = NaN;
    return;
end

switch statistics.family_id
    case 'gaussian'
        mu = statistics.parameters(1);
        dW = 3.0 .* mu .* lambda_n;

    case 'nongaussian'
        mu = statistics.parameters(1);
        N = statistics.parameters(2);
        beta = inv_langevin(lambda_n ./ sqrt(N));
        dW = mu .* sqrt(N) .* beta;

    otherwise
        error('Micro_GenStrain:UnsupportedStatistics', ...
            'Unsupported chain statistics: %s.', statistics.family);
end
end

function out = is_admissible_stretch(lambda_n, statistics)
out = isfinite(lambda_n) && isreal(lambda_n) && lambda_n > 0.0;

if out && strcmp(statistics.family_id, 'nongaussian')
    N = statistics.parameters(2);
    out = lambda_n < sqrt(N);
end
end

function M = direction_tensor(theta, phi)
n = [cos(theta); sin(theta) .* cos(phi); sin(theta) .* sin(phi)];
M = n * n';
end

function [statisticsParameters, chainParameters, strainParameters] = ...
    split_parameters_by_family(parameters, statisticsFamily, ...
    chainFamily, strainFamily)
statisticsParameterCount = statistics_parameter_count(statisticsFamily);
chainParameterCount = strain_parameter_count(chainFamily);
strainParameterCount = strain_parameter_count(strainFamily);
[statisticsParameters, chainParameters, strainParameters] = ...
    split_parameters(parameters, statisticsParameterCount, ...
    chainParameterCount, strainParameterCount);
end

function [statisticsParameters, chainParameters, strainParameters] = ...
    split_parameters(parameters, statisticsParameterCount, ...
    chainParameterCount, strainParameterCount)
statisticsStart = 1;
statisticsStop = statisticsParameterCount;
chainStart = statisticsStop + 1;
chainStop = chainStart + chainParameterCount - 1;
strainStart = chainStop + 1;
strainStop = strainStart + strainParameterCount - 1;

statisticsParameters = parameters(statisticsStart:statisticsStop);
chainParameters = parameters(chainStart:chainStop);
strainParameters = parameters(strainStart:strainStop);
end

function count = statistics_parameter_count(statisticsFamily)
statisticsId = normalize_statistics_family(statisticsFamily);

switch statisticsId
    case 'gaussian'
        count = 1;
    case 'nongaussian'
        count = 2;
    otherwise
        error('Micro_GenStrain:UnsupportedStatistics', ...
            'Unsupported chain statistics: %s.', char(statisticsFamily));
end
end

function statistics = chain_statistics(statisticsFamily, parameters)
statisticsId = normalize_statistics_family(statisticsFamily);
parameters = parameters(:).';

switch statisticsId
    case 'gaussian'
        require_parameter_count(parameters, 1, 'Gaussian statistics');
        mu = parameters(1);

        if mu < 0.0
            error('Micro_GenStrain:InvalidParameter', ...
                'mu must be non-negative.');
        end

        statistics.family_id = 'gaussian';
        statistics.family = 'Gaussian';
        statistics.family_label = 'Gaussian';
        statistics.parameter_names = {'mu'};
        statistics.parameters = parameters;

    case 'nongaussian'
        require_parameter_count(parameters, 2, 'NonGaussian statistics');
        mu = parameters(1);
        N = parameters(2);

        if mu < 0.0
            error('Micro_GenStrain:InvalidParameter', ...
                'mu must be non-negative.');
        end
        if N <= 1.0
            error('Micro_GenStrain:InvalidParameter', ...
                'N must be greater than 1.');
        end

        statistics.family_id = 'nongaussian';
        statistics.family = 'NonGaussian';
        statistics.family_label = 'Non-Gaussian';
        statistics.parameter_names = {'mu', 'N'};
        statistics.parameters = parameters;

    otherwise
        error('Micro_GenStrain:UnsupportedStatistics', ...
            'Unsupported chain statistics: %s.', char(statisticsFamily));
end
end

function count = strain_parameter_count(family)
familyId = normalize_family(family);

switch familyId
    case {'hencky', 'biot'}
        count = 0;
    case {'sh', 'bi', 'cz'}
        count = 1;
    case {'cr', 'dn'}
        count = 2;
    otherwise
        error('Micro_GenStrain:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
end
end

function statisticsId = normalize_statistics_family(statisticsFamily)
statisticsId = lower(char(statisticsFamily));
statisticsId = strrep(statisticsId, '-', '_');
statisticsId = strrep(statisticsId, ' ', '_');

switch statisticsId
    case {'gaussian', 'gauss', 'g'}
        statisticsId = 'gaussian';
    case {'nongaussian', 'non_gaussian', 'ng'}
        statisticsId = 'nongaussian';
end
end

function familyId = normalize_family(family)
familyId = lower(char(family));
familyId = strrep(familyId, '-', '_');
familyId = strrep(familyId, ' ', '_');

switch familyId
    case {'seth_hill', 'sethhill'}
        familyId = 'sh';
    case {'he', 'log', 'logarithmic'}
        familyId = 'hencky';
    case {'bazant_itskov', 'bazantitskov'}
        familyId = 'bi';
    case {'curnier_rakotomanana', 'curnierrakotomanana'}
        familyId = 'cr';
    case {'curnier_zysset', 'curnierzysset'}
        familyId = 'cz';
    case {'darijani_naghdabadi', 'darijaninaghdabadi'}
        familyId = 'dn';
end
end

function names = prefixed_parameter_names(prefix, parameterNames)
names = cell(1, length(parameterNames));
for ii = 1:length(parameterNames)
    names{ii} = sprintf('%s_%s', prefix, parameterNames{ii});
end
end

function require_parameter_count(parameters, expectedCount, name)
if length(parameters) ~= expectedCount
    error('Micro_GenStrain:InvalidParameterCount', ...
        '%s expects %d parameter(s), got %d.', ...
        name, expectedCount, length(parameters));
end
end

function out = P(parameters, statisticsFamily, chainFamily, strainFamily, F)
out = incompressible_constraint( ...
    F * S(parameters, statisticsFamily, chainFamily, strainFamily, F), F);
end
