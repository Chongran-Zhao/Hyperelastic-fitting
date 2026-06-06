function model = Hill_GenStrain(parametersOrStrainFamily, strainFamily)
% Hill generalized strain model.
%
% Supported generalized strains:
%   SH, Hencky, Biot, BI, CR, CZ, DN.
%
% W = mu * E:E, with E = sum_a E(lambda_bar_a) * M_a.
%
% Reference:
%   Liu, J., Zhao, C., & Guan, J. (2025).
%   Modeling finite viscoelasticity based on the Green-Naghdi kinematic
%   assumption and generalized strains.
%   Journal of the Mechanics and Physics of Solids, 106346.
%
% parameters = [mu, strain parameters...]

if nargin < 1
    parametersOrStrainFamily = [];
end
if nargin < 2
    strainFamily = [];
end

[parameters, strainFamily, lowerBounds, upperBounds] = ...
    parse_inputs(nargin, parametersOrStrainFamily, strainFamily);

strainParameterCount = strain_parameter_count(strainFamily);
expectedParameterCount = 1 + strainParameterCount;

if length(parameters) ~= expectedParameterCount
    error('Hill_GenStrain:InvalidParameters', ...
        'The parameters must be [mu, strain parameters...], expected %d entries.', ...
        expectedParameterCount);
end

mu = parameters(1);
strainParameters = parameters(2:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

if mu < 0.0
    error('Hill_GenStrain:InvalidParameter', 'mu must be non-negative.');
end

model.name = sprintf('Hill-%s', strainInfo.family);
model.parameters = parameters;
model.lower_bounds = lowerBounds;
model.upper_bounds = upperBounds;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu'}, strainInfo.parameter_names];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Hill_GenStrain(parameters, strainFamily);
end

function [parameters, strainFamily, lowerBounds, upperBounds] = ...
    parse_inputs(inputCount, parametersOrStrainFamily, strainFamily)
if inputCount < 1 || isempty(parametersOrStrainFamily)
    error('Hill_GenStrain:InvalidCall', ...
        ['Use Hill_GenStrain(strainFamily) for default parameters or ', ...
        'Hill_GenStrain(parameters, strainFamily) for explicit parameters.']);
end

if is_text(parametersOrStrainFamily)
    if inputCount ~= 1
        error('Hill_GenStrain:InvalidCall', ...
            'Use Hill_GenStrain(strainFamily) for default parameters.');
    end

    strainFamily = parametersOrStrainFamily;
    [parameters, lowerBounds, upperBounds] = default_parameters(strainFamily);
    return;
end

parameters = parametersOrStrainFamily(:).';
if inputCount < 2 || isempty(strainFamily)
    strainFamily = 'SH';
end
[~, lowerBounds, upperBounds] = default_parameters(strainFamily);
end

function out = is_text(value)
out = ischar(value) || ...
    ((exist('isstring', 'builtin') || exist('isstring', 'file')) && isstring(value));
end

function [parameters, lowerBounds, upperBounds] = default_parameters(strainFamily)
[strainParameters, strainLowerBounds, strainUpperBounds] = ...
    default_family_parameters(strainFamily);

parameters = [0.01, strainParameters];
lowerBounds = [0.0, strainLowerBounds];
upperBounds = [Inf, strainUpperBounds];
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
        error('Hill_GenStrain:UnsupportedFamily', ...
            'Unsupported generalized strain family: %s.', char(family));
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E = spectral_tensor(strain.values, kin.V);
W = mu .* contract(E, E);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);
E = spectral_tensor(strain.values, kin.V);
T = 2.0 .* mu .* E;
S_bar = contract(T, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, strainFamily, F)
out = incompressible_constraint(F * S(parameters, strainFamily, F), F);
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
        error('Hill_GenStrain:UnsupportedFamily', ...
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
