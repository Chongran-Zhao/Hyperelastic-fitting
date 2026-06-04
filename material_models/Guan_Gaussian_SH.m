function model = Guan_Gaussian_SH(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Seth-Hill chain strain.
%
% The chain scale function E_hat is the Seth-Hill strain:
%   E_hat(lambda_n) = (lambda_n^m_hat - 1)/m_hat.
% The limit m_hat = 0 gives the Hencky chain strain.
%
% The macroscopic generalized strain E_bar is selectable. For
% Gaussian-chain statistics, the sphere average is evaluated by Lebedev
% quadrature:
%   Psi = 3/2*mu<(((1 + m_hat*E_bar:(n*n'))^2)^(1/m_hat)) - 1>.
% This expression is used on the physical branch
%   1 + m_hat*E_bar:(n*n') > 0.
% The reference constant is chosen so that Psi = 0 at F = I.
%
% The isochoric second Piola-Kirchhoff stress follows from
%   S = J^(-2/3) DEV(S_bar),
%   S_bar = dPsi_dE_bar : Q_bar,
% where Q_bar = 2*dE_bar/dC_bar is Hill's fourth-order projection tensor.
%
% Supported macroscopic generalized strains:
%   SH, Hencky, Biot, CR, CZ, DN.
%
% Reference:
%   Guan, J., Li, X., Yuan, H., & Liu, J. (2025).
%   Hyperelastic modeling based on generalized Landau invariants and
%   multi-stage calibration.
%   Journal of the Mechanics and Physics of Solids, 106338.
%
% parameters = [mu, m_hat, strain parameters...]

if nargin < 2 || isempty(strainFamily)
    strainFamily = 'SH';
end

parameters = parameters(:).';
if length(parameters) < 2
    error('Guan_Gaussian_SH:InvalidParameters', ...
        'The parameters must be [mu, m_hat, strain parameters...].');
end

mu = parameters(1);
mHat = parameters(2);
strainParameters = parameters(3:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_SH', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.chain_strain_family = 'Seth-Hill';
model.chain_strain_parameter = mHat;
model.parameter_names = [{'mu', 'm_hat'}, ...
    prefix_strain_parameter_names(strainInfo.parameter_names)];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_SH(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_SH:InvalidParameter', 'mu must be non-negative.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
mHat = parameters(2);
strainParameters = parameters(3:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

energyDensity = @(theta, phi) chain_energy_density( ...
    directional_strain(E_bar, theta, phi), mHat);
W = 1.5 .* mu .* lebedev_quadrature(energyDensity);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
mHat = parameters(2);
strainParameters = parameters(3:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);
E_bar = spectral_tensor(strain.values, kin.V);

dPsi_dE_bar = 3.0 .* mu .* lebedev_quadrature( ...
    @(theta, phi) chain_energy_derivative( ...
    directional_strain(E_bar, theta, phi), mHat) .* ...
    direction_projection(theta, phi));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function value = chain_energy_density(directionalE, mHat)
if abs(mHat) < 1.0e-12
    value = exp(2.0 .* directionalE) - 1.0;
    return;
end

base = 1.0 + mHat .* directionalE;
if base <= 0.0
    value = NaN;
    return;
end

value = (base .^ 2.0) .^ (1.0 ./ mHat) - 1.0;
end

function value = chain_energy_derivative(directionalE, mHat)
if abs(mHat) < 1.0e-12
    value = exp(2.0 .* directionalE);
    return;
end

base = 1.0 + mHat .* directionalE;
if base <= 0.0
    value = NaN;
    return;
end

value = base .* (base .^ 2.0) .^ (1.0 ./ mHat - 1.0);
end

function value = directional_strain(E, theta, phi)
M = direction_projection(theta, phi);
value = contract(E, M);
end

function M = direction_projection(theta, phi)
n = [cos(theta); sin(theta) .* cos(phi); sin(theta) .* sin(phi)];
M = n * n';
end

function names = prefix_strain_parameter_names(names)
for ii = 1:length(names)
    names{ii} = ['strain_', names{ii}];
end
end

function out = P(parameters, strainFamily, F)
out = incompressible_constraint(F * S(parameters, strainFamily, F), F);
end
