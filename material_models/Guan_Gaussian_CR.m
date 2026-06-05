function model = Guan_Gaussian_CR(parameters, strainFamily)
% Guan Gaussian-chain micro-macro model with Curnier-Rakotomanana chain strain.
%
% The chain scale function E_hat is the Curnier-Rakotomanana strain:
%   E_hat(lambda_n) = (lambda_n^m_hat - lambda_n^(-n_hat))/(m_hat+n_hat).
% The chain-scale parameters satisfy m_hat*n_hat > 0.
%
% The macroscopic generalized strain E_bar is selectable. With
% xi = E_bar:(n*n'), lambda_n is implicitly determined by
%   lambda_n^m_hat - lambda_n^(-n_hat) = (m_hat+n_hat)*xi.
% Gaussian-chain statistics then give:
%   Psi = 3/2*mu*<lambda_n^2 - 1>.
% The reference constant is chosen so that Psi = 0 at F = I.
%
% The isochoric second Piola-Kirchhoff stress follows from
%   S = J^(-2/3) DEV(S_bar),
%   S_bar = dPsi_dE_bar : Q_bar,
% where Q_bar = 2*dE_bar/dC_bar is Hill's fourth-order projection tensor.
%
% Supported macroscopic generalized strains:
%   SH, Hencky, Biot, BI, CR, CZ, DN.
%
% Reference:
%   Guan, J., Li, X., Yuan, H., & Liu, J. (2025).
%   Hyperelastic modeling based on generalized Landau invariants and
%   multi-stage calibration.
%   Journal of the Mechanics and Physics of Solids, 106338.
%
% parameters = [mu, m_hat, n_hat, strain parameters...]

if nargin < 2 || isempty(strainFamily)
    strainFamily = 'SH';
end

parameters = parameters(:).';
if length(parameters) < 3
    error('Guan_Gaussian_CR:InvalidParameters', ...
        'The parameters must be [mu, m_hat, n_hat, strain parameters...].');
end

mu = parameters(1);
mHat = parameters(2);
nHat = parameters(3);
strainParameters = parameters(4:end);
strainInfo = generalized_strain(strainFamily, ones(3, 1), strainParameters);

model.name = sprintf('Guan-Gaussian-%s_CR', strainInfo.family);
model.parameters = parameters;
model.strain_family = strainInfo.family;
model.strain_family_label = strainInfo.family_label;
model.chain_strain_family = 'Curnier-Rakotomanana';
model.chain_strain_parameters = [mHat, nHat];
model.parameter_names = [{'mu', 'm_hat', 'n_hat'}, ...
    prefix_strain_parameter_names(strainInfo.parameter_names)];
model.energy = @(F) Energy(parameters, strainFamily, F);
model.S = @(F) S(parameters, strainFamily, F);
model.P = @(F) P(parameters, strainFamily, F);
model.set_parameters = @(parameters) Guan_Gaussian_CR(parameters, strainFamily);

if mu < 0.0
    error('Guan_Gaussian_CR:InvalidParameter', 'mu must be non-negative.');
end
if mHat .* nHat <= 0.0
    error('Guan_Gaussian_CR:InvalidParameter', ...
        'm_hat and n_hat must satisfy m_hat*n_hat > 0.');
end
end

function W = Energy(parameters, strainFamily, F)
mu = parameters(1);
mHat = parameters(2);
nHat = parameters(3);
strainParameters = parameters(4:end);

kin = kinematics(F, 'F');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters);
E_bar = spectral_tensor(strain.values, kin.V);

energyDensity = @(theta, phi) lambda_n_squared( ...
    contract(E_bar, direction_tensor(theta, phi)), mHat, nHat) - 1.0;
W = 1.5 .* mu .* lebedev_quadrature(energyDensity);
end

function out = S(parameters, strainFamily, F)
mu = parameters(1);
mHat = parameters(2);
nHat = parameters(3);
strainParameters = parameters(4:end);

C = F' * F;
kin = kinematics(C, 'C');
strain = generalized_strain(strainFamily, kin.lambda_bar, strainParameters, kin.V);
E_bar = spectral_tensor(strain.values, kin.V);

dPsi_dE_bar = 1.5 .* mu .* lebedev_quadrature( ...
    @(theta, phi) d_lambda_n_squared_d_xi( ...
    contract(E_bar, direction_tensor(theta, phi)), mHat, nHat) .* ...
    direction_tensor(theta, phi));
S_bar = contract(dPsi_dE_bar, strain.Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function value = lambda_n_squared(xi, mHat, nHat)
lambda = lambda_n(xi, mHat, nHat);
value = lambda .* lambda;
end

function value = d_lambda_n_squared_d_xi(xi, mHat, nHat)
lambda = lambda_n(xi, mHat, nHat);
denominator = mHat .* lambda .^ (mHat - 1.0) + ...
    nHat .* lambda .^ (-nHat - 1.0);
value = 2.0 .* (mHat + nHat) .* lambda ./ denominator;
end

function lambda = lambda_n(xi, mHat, nHat)
[lower, upper] = log_stretch_bracket(xi, mHat, nHat);
u = 0.5 .* (lower + upper);
gLower = cr_residual(lower, xi, mHat, nHat);

for iteration = 1:60
    [gValue, dgValue] = cr_residual_and_derivative(u, xi, mHat, nHat);
    if abs(gValue) <= 1.0e-12
        lambda = exp(u);
        return;
    end

    if same_sign(gValue, gLower)
        lower = u;
        gLower = gValue;
    else
        upper = u;
    end

    newtonStep = gValue ./ dgValue;
    uCandidate = u - newtonStep;
    if ~isfinite(uCandidate) || uCandidate <= lower || uCandidate >= upper
        uCandidate = 0.5 .* (lower + upper);
    end

    if abs(uCandidate - u) <= 1.0e-12
        lambda = exp(uCandidate);
        return;
    end
    u = uCandidate;
end

lambda = exp(u);
end

function [lower, upper] = log_stretch_bracket(xi, mHat, nHat)
lower = -1.0;
upper = 1.0;

for iteration = 1:80
    gLower = cr_residual(lower, xi, mHat, nHat);
    gUpper = cr_residual(upper, xi, mHat, nHat);
    if gLower == 0.0
        upper = lower;
        return;
    end
    if gUpper == 0.0
        lower = upper;
        return;
    end
    if ~same_sign(gLower, gUpper)
        return;
    end

    lower = 2.0 .* lower;
    upper = 2.0 .* upper;
end

error('Guan_Gaussian_CR:RootNotBracketed', ...
    'Unable to bracket lambda_n for xi = %.12g.', xi);
end

function [value, derivative] = cr_residual_and_derivative(u, xi, mHat, nHat)
expM = exp(mHat .* u);
expN = exp(-nHat .* u);
value = expM - expN - (mHat + nHat) .* xi;
derivative = mHat .* expM + nHat .* expN;
end

function value = cr_residual(u, xi, mHat, nHat)
value = cr_residual_and_derivative(u, xi, mHat, nHat);
end

function out = same_sign(a, b)
out = (a < 0.0 && b < 0.0) || (a > 0.0 && b > 0.0);
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
