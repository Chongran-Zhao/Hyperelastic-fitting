function model = Zhao_NonGaussian(parametersOrChainFamily, chainFamilyOrStrainFamily, strainFamily)
% General Zhao non-Gaussian chain-network model.
%
% The chain strain E_hat(lambda_n) and the macroscopic generalized strain
% E_bar are independently selectable. For each quadrature direction n,
%   E_hat(lambda_n) = E_bar:(n*n')
% and lambda_n is obtained from the inverse of E_hat.
%
% The isochoric second Piola-Kirchhoff stress follows from
%   S = J^(-2/3) DEV(S_bar),
%   S_bar = dPsi/dE_bar : Q_bar,
% where Q_bar = 2*dE_bar/dC_bar is Hill's fourth-order projection tensor.
%
% Supported chain and macroscopic generalized strains:
%   SH, Hencky, Biot, BI, CR, CZ, DN.
%
% Examples:
%   Zhao_NonGaussian('Hencky', 'Hencky')
%   Zhao_NonGaussian([1.0, 50.0, 1.0], 'BI', 'Biot')
%
% Explicit parameter vector:
%   parameters = [mu, N, chain parameters..., strain parameters...]

if nargin < 1
    parametersOrChainFamily = [];
end
if nargin < 2
    chainFamilyOrStrainFamily = [];
end
if nargin < 3
    strainFamily = [];
end

[parameters, chainFamily, strainFamily, lowerBounds, upperBounds] = ...
    parse_inputs(nargin, parametersOrChainFamily, ...
    chainFamilyOrStrainFamily, strainFamily);
chainParameterCount = strain_parameter_count(chainFamily);
strainParameterCount = strain_parameter_count(strainFamily);
expectedParameterCount = 2 + chainParameterCount + strainParameterCount;

if length(parameters) ~= expectedParameterCount
    error('Zhao_NonGaussian:InvalidParameters', ...
        ['The parameters must be [mu, N, chain parameters..., ', ...
        'strain parameters...], expected %d entries.'], expectedParameterCount);
end

mu = parameters(1);
N = parameters(2);
[chainParameters, strainParameters] = split_parameters( ...
    parameters, chainParameterCount, strainParameterCount);

chainInfo = generalized_strain(chainFamily, ones(3, 1), chainParameters);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Zhao-NonGaussian-Ehat%s-E%s', ...
    chainInfo.family, strainInfo.family);
model.parameters = parameters;
model.lower_bounds = lowerBounds;
model.upper_bounds = upperBounds;
model.chain_family = chainInfo.family;
model.chain_family_label = chainInfo.family_label;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu', 'N'}, ...
    prefixed_parameter_names('chain', chainInfo.parameter_names), ...
    prefixed_parameter_names('strain', strainInfo.parameter_names)];
model.energy = @(F) Energy(parameters, chainFamily, strainFamily, F);
model.S = @(F) S(parameters, chainFamily, strainFamily, F);
model.P = @(F) P(parameters, chainFamily, strainFamily, F);
model.set_parameters = @(parameters) Zhao_NonGaussian( ...
    parameters, chainFamily, strainFamily);

if mu < 0.0
    error('Zhao_NonGaussian:InvalidParameter', 'mu must be non-negative.');
end
if N <= 1.0
    error('Zhao_NonGaussian:InvalidParameter', 'N must be greater than 1.');
end
end

function [parameters, chainFamily, strainFamily, lowerBounds, upperBounds] = ...
    parse_inputs(inputCount, parametersOrChainFamily, ...
    chainFamilyOrStrainFamily, strainFamily)
if inputCount < 1 || isempty(parametersOrChainFamily)
    chainFamily = 'Hencky';
    strainFamily = 'Hencky';
    [parameters, lowerBounds, upperBounds] = default_parameters( ...
        chainFamily, strainFamily);
    return;
end

if is_text(parametersOrChainFamily)
    chainFamily = parametersOrChainFamily;
    if inputCount < 2 || isempty(chainFamilyOrStrainFamily)
        strainFamily = 'Hencky';
    else
        strainFamily = chainFamilyOrStrainFamily;
    end
    [parameters, lowerBounds, upperBounds] = default_parameters( ...
        chainFamily, strainFamily);
    return;
end

parameters = parametersOrChainFamily(:).';
if inputCount < 2 || isempty(chainFamilyOrStrainFamily)
    chainFamily = 'Hencky';
else
    chainFamily = chainFamilyOrStrainFamily;
end
if inputCount < 3 || isempty(strainFamily)
    strainFamily = 'Hencky';
end
[~, lowerBounds, upperBounds] = default_parameters(chainFamily, strainFamily);
end

function out = is_text(value)
out = ischar(value) || ...
    ((exist('isstring', 'builtin') || exist('isstring', 'file')) && isstring(value));
end

function [parameters, lowerBounds, upperBounds] = default_parameters( ...
    chainFamily, strainFamily)
[chainParameters, chainLowerBounds, chainUpperBounds] = ...
    default_family_parameters(chainFamily);
[strainParameters, strainLowerBounds, strainUpperBounds] = ...
    default_family_parameters(strainFamily);

parameters = [1.0, 50.0, chainParameters, strainParameters];
lowerBounds = [0.0, 1.0 + 1.0e-8, chainLowerBounds, strainLowerBounds];
upperBounds = [Inf, Inf, chainUpperBounds, strainUpperBounds];
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
        error('Zhao_NonGaussian:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
end
end

function W = Energy(parameters, chainFamily, strainFamily, F)
mu = parameters(1);
N = parameters(2);
[chainParameters, strainParameters] = split_parameters_by_family( ...
    parameters, chainFamily, strainFamily);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
chain = generalized_strain(chainFamily, ones(3, 1), chainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

W0 = chain_energy_density(1.0, N);
energyDensity = @(theta, phi) ...
    chain_energy_density(chain_stretch(E_bar, chain, theta, phi), N) - W0;
W = mu .* N .* lebedev_quadrature(energyDensity);
end

function out = S(parameters, chainFamily, strainFamily, F)
mu = parameters(1);
N = parameters(2);
[chainParameters, strainParameters] = split_parameters_by_family( ...
    parameters, chainFamily, strainFamily);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);
chain = generalized_strain(chainFamily, ones(3, 1), chainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

dPsi_dE_bar = mu .* sqrt(N) .* lebedev_quadrature( ...
    @(theta, phi) chain_stress_weight( ...
    E_bar, chain, chainFamily, chainParameters, theta, phi, N) .* ...
    direction_tensor(theta, phi));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function lambda_n = chain_stretch(E_bar, chain, theta, phi)
xi_n = contract(E_bar, direction_tensor(theta, phi));
lambda_n = chain.inverse(xi_n);
lambda_n = lambda_n(1);

if ~isreal(lambda_n) || lambda_n <= 0.0
    lambda_n = NaN;
end
end

function value = chain_stress_weight(E_bar, chain, chainFamily, ...
    chainParameters, theta, phi, N)
lambda_n = chain_stretch(E_bar, chain, theta, phi);
if ~isfinite(lambda_n) || lambda_n >= sqrt(N)
    value = NaN;
    return;
end

beta = inv_langevin(lambda_n ./ sqrt(N));
chainAtStretch = generalized_strain(chainFamily, lambda_n, chainParameters);
value = beta ./ chainAtStretch.derivatives(1);
end

function W = chain_energy_density(lambda_n, N)
if ~isfinite(lambda_n) || lambda_n <= 0.0 || lambda_n >= sqrt(N)
    W = NaN;
    return;
end

beta = inv_langevin(lambda_n ./ sqrt(N));
W = lambda_n .* beta ./ sqrt(N) + log(beta ./ sinh(beta));
end

function M = direction_tensor(theta, phi)
n = [cos(theta); sin(theta) .* cos(phi); sin(theta) .* sin(phi)];
M = n * n';
end

function [chainParameters, strainParameters] = split_parameters_by_family( ...
    parameters, chainFamily, strainFamily)
chainParameterCount = strain_parameter_count(chainFamily);
strainParameterCount = strain_parameter_count(strainFamily);
[chainParameters, strainParameters] = split_parameters( ...
    parameters, chainParameterCount, strainParameterCount);
end

function [chainParameters, strainParameters] = split_parameters( ...
    parameters, chainParameterCount, strainParameterCount)
chainStart = 3;
chainStop = chainStart + chainParameterCount - 1;
strainStart = chainStop + 1;
strainStop = strainStart + strainParameterCount - 1;

chainParameters = parameters(chainStart:chainStop);
strainParameters = parameters(strainStart:strainStop);
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
        error('Zhao_NonGaussian:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
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

function out = P(parameters, chainFamily, strainFamily, F)
out = incompressible_constraint( ...
    F * S(parameters, chainFamily, strainFamily, F), F);
end
