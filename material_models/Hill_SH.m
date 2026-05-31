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
kin = kinematics(C, 'C');
strain = Strain_values(parameters, kin.lambda_bar);
E = spectral_tensor(strain.values, kin.V);
W = mu .* contract(E, E);
end

function out = S(parameters, F)
mu = parameters(1);

C = F' * F;
kin = kinematics(C, 'C');
strain = Strain_values(parameters, kin.lambda_bar);
E = spectral_tensor(strain.values, kin.V);
T = 2.0 .* mu .* E;
Q = hill_Q_proj(kin.lambda_bar, strain.values, strain.derivatives, kin.V);
S_bar = contract(T, Q);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end

function out = Strain_values(parameters, lambda)
m = parameters(2);
out.values = (lambda .^ m - 1.0) ./ m;
out.derivatives = lambda .^ (m - 1.0);
end
