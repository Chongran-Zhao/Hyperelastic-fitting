function model = Neo_Hookean(parameters)
% Incompressible Neo-Hookean model.
% I1_bar = sum(lambda_bar.^2).
% W = 0.5*mu*(I1_bar - 3).
% parameters = [mu]
parameters = parameters(:).';

model.name = 'Neo-Hookean';
model.parameters = parameters;
model.parameter_names = {'mu'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Neo_Hookean(parameters);
end

function W = Energy(parameters, F)
mu = parameters(1);

kinematics = Kinematics(F, 'F');
I1 = sum(kinematics.lambda_bar .^ 2);
W = 0.5 .* mu .* (I1 - 3.0);
end

function out = S(parameters, F)
mu = parameters(1);

C = F' * F;
kinematics = Kinematics(C, 'C');
S_bar = mu .* eye(3);
out = kinematics.J^(-2.0 / 3.0) .* Dev(S_bar, C);
end

function out = P(parameters, F)
out = Incompressible_constraint(F * S(parameters, F), F);
end
