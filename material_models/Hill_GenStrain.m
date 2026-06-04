function model = Hill_GenStrain(parameters, strainFamily)
% Hill generalized strain model.
%
% Supported generalized strains:
%   SH, Hencky, Biot, CR, CZ, DN.
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

if nargin < 2 || isempty(strainFamily)
    strainFamily = 'SH';
end

parameters = parameters(:).';
if isempty(parameters)
    error('Hill_GenStrain:InvalidParameters', ...
        'At least the modulus parameter mu is required.');
end

mu = parameters(1);
strainParameters = parameters(2:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

if mu < 0.0
    error('Hill_GenStrain:InvalidParameter', 'mu must be non-negative.');
end

model.name = sprintf('Hill-%s', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu'}, strainInfo.parameter_names];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Hill_GenStrain(parameters, strainFamily);
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
