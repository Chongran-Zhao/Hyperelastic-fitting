function model = Zhan_Gaussian(parameters)
% Zhan Gaussian chain-network model.
% W = mu*(sum(lambda_bar)^2 + 2*sum(lambda_bar.^2)).
% parameters = [mu]
parameters = parameters(:).';

model.name = 'Zhan-Gaussian';
model.parameters = parameters;
model.parameter_names = {'mu'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Zhan_Gaussian(parameters);
end

function W = Energy(parameters, F)
kinematics = Kinematics(F, 'F');
mu = parameters(1);
lambda = kinematics.lambda_bar;
W = mu .* (sum(lambda) .^ 2 + 2.0 .* sum(lambda .^ 2));
end

function out = S(parameters, F)
C = F' * F;
kinematics = Kinematics(C, 'C');
mu = parameters(1);
lambda = kinematics.lambda_bar;

values = 2.0 .* mu .* ((sum(lambda) ./ lambda) + 2.0);
S_bar = Spectral_tensor(values, kinematics.V);
out = kinematics.J^(-2.0 / 3.0) .* Dev(S_bar, C);
end

function out = P(parameters, F)
out = Incompressible_constraint(F * S(parameters, F), F);
end
