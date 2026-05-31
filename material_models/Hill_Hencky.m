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
kin = kinematics(C, 'C');
strain = Strain_values(kin.lambda_bar);
E = spectral_tensor(strain.values, kin.V);
W = mu .* contract(E, E);
end

function out = S(parameters, F)
mu = parameters(1);

C = F' * F;
kin = kinematics(C, 'C');
strain = Strain_values(kin.lambda_bar);
E = spectral_tensor(strain.values, kin.V);
T = 2.0 .* mu .* E;
Q = hill_Q_proj(kin.lambda_bar, strain.values, strain.derivatives, kin.V);
S_bar = contract(T, Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end

function out = Strain_values(lambda)
out.values = log(lambda);
out.derivatives = 1.0 ./ lambda;
end
