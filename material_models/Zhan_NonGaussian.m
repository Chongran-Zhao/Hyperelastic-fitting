function model = Zhan_NonGaussian(parameters)
% Zhan non-Gaussian chain-network model.
% Directional chain stretch is integrated over the unit sphere.
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
lambda = kin.lambda_bar;

lambda_i = @(theta, phi) lambda(1) .* cos(theta) .^ 2 + ...
    lambda(2) .* sin(theta) .^ 2 .* cos(phi) .^ 2 + ...
    lambda(3) .* sin(theta) .^ 2 .* sin(phi) .^ 2;
beta = @(theta, phi) inv_langevin(lambda_i(theta, phi) ./ sqrt(N));

Psi = @(theta, phi) mu .* N .* ...
    (lambda_i(theta, phi) .* beta(theta, phi) ./ sqrt(N) + ...
    log(beta(theta, phi) ./ sinh(beta(theta, phi)))) .* sin(theta);
W = integral2(Psi, 0, pi, 0, 2.0 .* pi, 'AbsTol', 1.0e-12, 'RelTol', 1.0e-9);
end

function out = S(parameters, F)
C = F' * F;
kin = kinematics(C, 'C');
mu = parameters(1);
N = parameters(2);
lambda = kin.lambda_bar;

lambda_i = @(theta, phi) lambda(1) .* cos(theta) .^ 2 + ...
    lambda(2) .* sin(theta) .^ 2 .* cos(phi) .^ 2 + ...
    lambda(3) .* sin(theta) .^ 2 .* sin(phi) .^ 2;
beta = @(theta, phi) inv_langevin(lambda_i(theta, phi) ./ sqrt(N));

temp_1 = @(theta, phi) beta(theta, phi) .* cos(theta) .* cos(theta);
temp_2 = @(theta, phi) beta(theta, phi) .* sin(theta) .* sin(theta) .* ...
    cos(phi) .* cos(phi);
temp_3 = @(theta, phi) beta(theta, phi) .* sin(theta) .* sin(theta) .* ...
    sin(phi) .* sin(phi);

values = zeros(3, 1);
values(1) = mu .* sqrt(N) ./ lambda(1) .* lebedev_quadrature(temp_1);
values(2) = mu .* sqrt(N) ./ lambda(2) .* lebedev_quadrature(temp_2);
values(3) = mu .* sqrt(N) ./ lambda(3) .* lebedev_quadrature(temp_3);

S_bar = spectral_tensor(values, kin.V);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
