function model = Hill_Hencky(parameters)
% Hill generalized strain model with Hencky logarithmic strain.
% E_i = log(lambda_bar_i).
% parameters = [mu]
parameters = parameters(:).';

model.name = 'Hill-Hencky';
model.parameters = parameters;
model.parameter_names = {'mu'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Hill_Hencky(parameters);
end

function W = Energy(parameters, F)
mu = parameters(1);

C = F' * F;
kinematics = Kinematics(C, 'C');
strain = Strain_values(kinematics.lambda_bar);
E = Spectral_tensor(strain.values, kinematics.V);
W = mu .* Contract(E, E);
end

function out = S(parameters, F)
mu = parameters(1);

C = F' * F;
kinematics = Kinematics(C, 'C');
strain = Strain_values(kinematics.lambda_bar);
E = Spectral_tensor(strain.values, kinematics.V);
T = 2.0 .* mu .* E;
Q = Hill_Q_proj(kinematics.lambda_bar, strain.values, strain.derivatives, kinematics.V);
S_bar = Contract(T, Q);
out = kinematics.J^(-2.0 / 3.0) .* Dev(S_bar, C);
end

function out = P(parameters, F)
out = Incompressible_constraint(F * S(parameters, F), F);
end

function out = Strain_values(lambda)
out.values = log(lambda);
out.derivatives = 1.0 ./ lambda;
end
