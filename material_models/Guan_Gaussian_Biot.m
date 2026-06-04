function model = Guan_Gaussian_Biot(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Biot chain strain.
%
% The chain scale function E_hat is fixed to the Biot strain:
%   E_hat(lambda_n) = lambda_n - 1.
%
% The macroscopic generalized strain E_bar is selectable. For
% Gaussian-chain statistics, the sphere average gives:
%   Psi = mu/10*((tr(E_bar))^2 + 2*E_bar:E_bar + 10*tr(E_bar)).
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
    error('Guan_Gaussian_Biot:InvalidParameters', ...
        'At least the modulus parameter mu is required.');
end

mu = parameters(1);
strainParameters = parameters(2:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_Biot', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu'}, strainInfo.parameter_names];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_Biot(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_Biot:InvalidParameter', 'mu must be non-negative.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);
traceE = trace(E_bar);
W = mu ./ 10.0 .* (traceE .* traceE + 2.0 .* contract(E_bar, E_bar) + ...
    10.0 .* traceE);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);

E_bar = spectral_tensor(strain.values, kin.V);
dPsi_dE_bar = mu ./ 10.0 .* ...
    (2.0 .* trace(E_bar) .* eye(size(E_bar)) + 4.0 .* E_bar + 10.0 .* eye(size(E_bar)));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, strainFamily, F)
out = incompressible_constraint(F * S(parameters, strainFamily, F), F);
end
