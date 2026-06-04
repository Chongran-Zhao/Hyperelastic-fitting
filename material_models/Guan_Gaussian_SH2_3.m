function model = Guan_Gaussian_SH2_3(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Seth-Hill 2/3 chain strain.
%
% The chain scale function E_hat is fixed to the Seth-Hill strain with
% m = 2/3:
%   E_hat(lambda_n) = 3/2*(lambda_n^(2/3) - 1).
%
% The macroscopic generalized strain E_bar is selectable. For
% Gaussian-chain statistics, the sphere average gives:
%   Psi = mu/945*(4*I1^3 + 126*I1^2 + 810*I1 + 24*I1*I2
%                 + 252*I2 + 32*I3),
% where I1 = tr(E_bar), I2 = tr(E_bar^2), and I3 = tr(E_bar^3).
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
    error('Guan_Gaussian_SH2_3:InvalidParameters', ...
        'At least the modulus parameter mu is required.');
end

mu = parameters(1);
strainParameters = parameters(2:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_SH2_3', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.parameter_names = [{'mu'}, strainInfo.parameter_names];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_SH2_3(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_SH2_3:InvalidParameter', 'mu must be non-negative.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);
[I1, I2, I3] = landau_invariants(E_bar);

W = mu ./ 945.0 .* (4.0 .* I1 .^ 3.0 + 126.0 .* I1 .^ 2.0 + ...
    810.0 .* I1 + 24.0 .* I1 .* I2 + 252.0 .* I2 + 32.0 .* I3);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
strainParameters = parameters(2:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);

E_bar = spectral_tensor(strain.values, kin.V);
[I1, I2, ~] = landau_invariants(E_bar);
I = eye(size(E_bar));
dPsi_dE_bar = mu ./ 945.0 .* ...
    ((12.0 .* I1 .^ 2.0 + 252.0 .* I1 + 810.0 + 24.0 .* I2) .* I + ...
    (48.0 .* I1 + 504.0) .* E_bar + 96.0 .* (E_bar * E_bar));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function [I1, I2, I3] = landau_invariants(E)
I1 = trace(E);
I2 = contract(E, E);
I3 = trace(E * E * E);
end

function out = P(parameters, strainFamily, F)
out = incompressible_constraint(F * S(parameters, strainFamily, F), F);
end
