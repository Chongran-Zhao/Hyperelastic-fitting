function model = Guan_Gaussian_CZ(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Curnier-Zysset chain strain.
%
% The chain scale function E_hat is the Curnier-Zysset strain:
%   E_hat(lambda_n) = (lambda_n^m_hat - lambda_n^(-m_hat))/(2*m_hat).
% The chain-scale parameter m_hat is nonzero.
%
% The macroscopic generalized strain E_bar is selectable. With
% xi = E_bar:(n*n'), Gaussian-chain statistics give:
%   Psi = 3/2*mu*<lambda_n^2 - 1>,
%   lambda_n^2 = exp(2/m_hat*asinh(m_hat*xi)).
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
    error('Guan_Gaussian_CZ:InvalidParameters', ...
        'The parameters must be [mu, m_hat, strain parameters...].');
end

mu = parameters(1);
mHat = parameters(2);
strainParameters = parameters(3:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_CZ', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.chain_strain_family = 'Curnier-Zysset';
model.chain_strain_parameter = mHat;
model.parameter_names = [{'mu', 'm_hat'}, ...
    prefix_strain_parameter_names(strainInfo.parameter_names)];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_CZ(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_CZ:InvalidParameter', 'mu must be non-negative.');
end
if abs(mHat) < 1.0e-12
    error('Guan_Gaussian_CZ:InvalidParameter', 'm_hat must be nonzero.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
mHat = parameters(2);
strainParameters = parameters(3:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

energyDensity = @(theta, phi) lambda_n_squared( ...
    contract(E_bar, direction_tensor(theta, phi)), mHat) - 1.0;
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

dPsi_dE_bar = 1.5 .* mu .* lebedev_quadrature( ...
    @(theta, phi) d_lambda_n_squared_d_xi( ...
    contract(E_bar, direction_tensor(theta, phi)), mHat) .* ...
    direction_tensor(theta, phi));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function value = lambda_n_squared(xi, mHat)
value = exp(2.0 ./ mHat .* asinh(mHat .* xi));
end

function value = d_lambda_n_squared_d_xi(xi, mHat)
argument = mHat .* xi;
value = 2.0 .* lambda_n_squared(xi, mHat) ./ ...
    sqrt(1.0 + argument .^ 2.0);
end

function M = direction_tensor(theta, phi)
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
