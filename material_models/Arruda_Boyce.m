function model = Arruda_Boyce(parameters)
% Arruda-Boyce eight-chain model with an inverse Langevin approximation.
% I1_bar = sum(lambda_bar.^2), and W is integrated from I1_bar = 3.
% parameters = [mu, N]
parameters = parameters(:).';

model.name = 'Arruda-Boyce';
model.parameters = parameters;
model.parameter_names = {'mu', 'N'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Arruda_Boyce(parameters);
end

function W = Energy(parameters, F)
kin = kinematics(F, 'F');
mu = parameters(1);
N = parameters(2);
I1 = sum(kin.lambda_bar .^ 2);

if abs(I1 - 3.0) < 1.0e-12
    W = 0.0;
    return;
end

integrand = @(x) 0.5 .* sqrt(N ./ (3.0 .* x)) .* mu .* ...
    inv_langevin(sqrt(x ./ (3.0 .* N)));
W = integral(integrand, 3.0, I1, 'AbsTol', 1.0e-10, 'RelTol', 1.0e-8);
end

function out = S(parameters, F)
C = F' * F;
kin = kinematics(C, 'C');
mu = parameters(1);
N = parameters(2);
I1 = sum(kin.lambda_bar .^ 2);
lambda_r = sqrt(I1 ./ (3.0 .* N));

value = sqrt(N ./ (3.0 .* I1)) .* mu .* inv_langevin(lambda_r);
S_bar = value .* eye(3);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
