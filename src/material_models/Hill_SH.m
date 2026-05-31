function model = Hill_SH(parameters)
% Hill generalized strain model with Seth-Hill strain.
% E_i = (lambda_bar_i^m - 1)/m, with m ~= 0.
% Use Hill_Hencky for the m = 0 logarithmic limit.
% parameters = [mu, m]
parameters = parameters(:).';

if abs(parameters(2)) < 1.0e-12
    error('Hill_SH:InvalidParameter', ...
        'Hill-SH requires m ~= 0. Use Hill_Hencky for the m = 0 limit.');
end

model.name = 'Hill-SH';
model.parameters = parameters;
model.parameter_names = {'mu', 'm'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Hill_SH(parameters);
end

function W = Energy(parameters, F)
mu = parameters(1);

C = F' * F;
kinematics = Kinematics(C, 'C');
strain = Strain_values(parameters, kinematics.lambda_bar);
E = Spectral_tensor(strain.values, kinematics.V);
W = mu .* Contract(E, E);
end

function out = S(parameters, F)
mu = parameters(1);

C = F' * F;
kinematics = Kinematics(C, 'C');
strain = Strain_values(parameters, kinematics.lambda_bar);
E = Spectral_tensor(strain.values, kinematics.V);
T = 2.0 .* mu .* E;
Q = Hill_Q_proj(kinematics.lambda_bar, strain.values, strain.derivatives, kinematics.V);
S_bar = Contract(T, Q);
out = kinematics.J^(-2.0 / 3.0) .* Dev(S_bar, C);
end

function out = P(parameters, F)
out = Incompressible_constraint(F * S(parameters, F), F);
end

function out = Strain_values(parameters, lambda)
m = parameters(2);
out.values = (lambda .^ m - 1.0) ./ m;
out.derivatives = lambda .^ (m - 1.0);
end
