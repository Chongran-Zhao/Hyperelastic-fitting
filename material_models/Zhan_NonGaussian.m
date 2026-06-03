function model = Zhan_NonGaussian(parameters)
% Zhan non-Gaussian chain-network model.
% Directional chain stretch is integrated over the unit sphere.
%
% Reference:
%   Zhan, L., Wang, S., Qu, S., Steinmann, P., & Xiao, R. (2023).
%   A new micro-macro transition for hyperelastic materials.
%   Journal of the Mechanics and Physics of Solids, 171, 105156.
%
% parameters = [mu, N]
parameters = parameters(:).';

model.name = 'Zhan-NonGaussian';
model.parameters = parameters;
model.parameter_names = {'mu', 'N'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Zhan_NonGaussian(parameters);
end

function W = Energy(parameters, F)
kin = kinematics(F, 'F');
mu = parameters(1);
N = parameters(2);
lambda_bar = kin.lambda_bar;

W0 = chain_energy_density(1.0, N);
energyDensity = @(theta, phi) ...
    chain_energy_density(chain_stretch(lambda_bar, theta, phi), N) - W0;
W = mu .* N .* lebedev_quadrature(energyDensity);
end

function out = S(parameters, F)
C = F' * F;
kin = kinematics(C, 'C');
mu = parameters(1);
N = parameters(2);
lambda_bar = kin.lambda_bar;

lambda_n = @(theta, phi) chain_stretch(lambda_bar, theta, phi);
beta = @(theta, phi) inv_langevin(lambda_n(theta, phi) ./ sqrt(N));

temp_1 = @(theta, phi) beta(theta, phi) .* cos(theta) .* cos(theta);
temp_2 = @(theta, phi) beta(theta, phi) .* sin(theta) .* sin(theta) .* ...
    cos(phi) .* cos(phi);
temp_3 = @(theta, phi) beta(theta, phi) .* sin(theta) .* sin(theta) .* ...
    sin(phi) .* sin(phi);

values = zeros(3, 1);
values(1) = mu .* sqrt(N) ./ lambda_bar(1) .* lebedev_quadrature(temp_1);
values(2) = mu .* sqrt(N) ./ lambda_bar(2) .* lebedev_quadrature(temp_2);
values(3) = mu .* sqrt(N) ./ lambda_bar(3) .* lebedev_quadrature(temp_3);

S_bar = spectral_tensor(values, kin.V);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function lambda_n = chain_stretch(lambda_bar, theta, phi)
lambda_n = lambda_bar(1) .* cos(theta) .^ 2 + ...
    lambda_bar(2) .* sin(theta) .^ 2 .* cos(phi) .^ 2 + ...
    lambda_bar(3) .* sin(theta) .^ 2 .* sin(phi) .^ 2;
end

function W = chain_energy_density(lambda_n, N)
beta = inv_langevin(lambda_n ./ sqrt(N));
W = lambda_n .* beta ./ sqrt(N) + log(beta ./ sinh(beta));
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
