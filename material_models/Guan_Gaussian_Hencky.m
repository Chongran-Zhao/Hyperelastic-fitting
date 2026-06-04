function model = Guan_Gaussian_Hencky(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Hencky chain strain.
%
% The chain scale function E_hat is fixed to the Hencky strain:
%   E_hat(lambda_n) = log(lambda_n).
%
% The macroscopic generalized strain E_bar is selectable. For
% Gaussian-chain statistics, the sphere average is evaluated by Lebedev
% quadrature:
%   Psi = 3/2*mu*<exp(2*E_bar:(n*n')) - 1>.
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
% parameters = [mu, strain parameters...]

if nargin < 2 || isempty(strainFamily)
    strainFamily = 'SH';
end

parameters = parameters(:).';
if isempty(parameters)
    error('Guan_Gaussian_Hencky:InvalidParameters', ...
        'At least the modulus parameter mu is required.');
end

mu = parameters(1);
strainParameters = parameters(2:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_Hencky', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu'}, strainInfo.parameter_names];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_Hencky(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_Hencky:InvalidParameter', 'mu must be non-negative.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

energyDensity = @(theta, phi) exp(2.0 .* directional_strain(E_bar, theta, phi)) - 1.0;
W = 1.5 .* mu .* lebedev_quadrature(energyDensity);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);
E_bar = spectral_tensor(strain.values, kin.V);

dPsi_dE_bar = 3.0 .* mu .* lebedev_quadrature( ...
    @(theta, phi) exp(2.0 .* directional_strain(E_bar, theta, phi)) .* ...
    direction_projection(theta, phi));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function value = directional_strain(E, theta, phi)
M = direction_projection(theta, phi);
value = contract(E, M);
end

function M = direction_projection(theta, phi)
n = [cos(theta); sin(theta) .* cos(phi); sin(theta) .* sin(phi)];
M = n * n';
end

function out = P(parameters, strainFamily, F)
out = incompressible_constraint(F * S(parameters, strainFamily, F), F);
end
