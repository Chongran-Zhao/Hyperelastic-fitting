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
out.values = (2.0 + m) ./ 8.0 .* lambda .^ 2.0 - ...
    (2.0 - m) ./ 8.0 .* lambda .^ (-2.0) - m ./ 4.0;
out.derivatives = (2.0 + m) ./ 4.0 .* lambda + ...
    (2.0 - m) ./ 4.0 .* lambda .^ (-3.0);
end
