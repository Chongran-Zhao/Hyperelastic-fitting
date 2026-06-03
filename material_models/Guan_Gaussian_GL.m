function model = Guan_Gaussian_GL(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Green-Lagrange chain strain.
%
% The chain scale function E_hat is fixed to the Green-Lagrange strain:
%   E_hat(lambda_n) = 1/2*(lambda_n^2 - 1).
%
% The macroscopic generalized strain E_bar is selectable. For
% Gaussian-chain statistics, the sphere average gives:
%   Psi = mu * tr(E_bar).
%
% The isochoric second Piola-Kirchhoff stress follows from
%   S = J^(-2/3) DEV(S_bar),
%   S_bar = dPsi_dE_bar : Q_bar,
% where Q_bar = 2*dE_bar/dC_bar is Hill's fourth-order projection tensor.
%
% Supported macroscopic generalized strains:
%   SH, Hencky, CR, CZ, DN.
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
    error('Guan_Gaussian_GL:InvalidParameters', ...
        'At least the modulus parameter mu is required.');
end

mu = parameters(1);
strainParameters = parameters(2:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_GL', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu'}, strainInfo.parameter_names];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_GL(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_GL:InvalidParameter', 'mu must be non-negative.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);
W = mu .* trace(E_bar);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);

E_bar = spectral_tensor(strain.values, kin.V);
dPsi_dE_bar = mu .* eye(size(E_bar));
Q_bar = hill_Q_proj(kin.lambda_bar, strain.values, strain.derivatives, kin.V);
S_bar = contract(dPsi_dE_bar, Q_bar);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, strainFamily, F)
out = incompressible_constraint(F * S(parameters, strainFamily, F), F);
end
