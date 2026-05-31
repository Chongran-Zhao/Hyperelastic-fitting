function model = Hill_CZ(parameters)
% Hill generalized strain model with Curnier-Zysset strain.
% E_i = (2 + m)/8*lambda_bar_i^2 - (2 - m)/8*lambda_bar_i^(-2) - m/4.
% parameters = [mu, m]
parameters = parameters(:).';

m = parameters(2);
if m < -2.0 || m > 2.0
    error('Hill_CZ:InvalidParameter', 'Hill-CZ requires -2 <= m <= 2.');
end

model.name = 'Hill-CZ';
model.parameters = parameters;
model.parameter_names = {'mu', 'm'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Hill_CZ(parameters);
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
out.values = (2.0 + m) ./ 8.0 .* lambda .^ 2.0 - ...
    (2.0 - m) ./ 8.0 .* lambda .^ (-2.0) - m ./ 4.0;
out.derivatives = (2.0 + m) ./ 4.0 .* lambda + ...
    (2.0 - m) ./ 4.0 .* lambda .^ (-3.0);
end
