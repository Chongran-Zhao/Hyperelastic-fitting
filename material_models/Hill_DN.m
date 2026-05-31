function model = Hill_DN(parameters)
% Hill generalized strain model with Darijani-Naghdabadi strain.
% E_i = (exp(m*(lambda_bar_i - 1)) - exp(n*(1/lambda_bar_i - 1)))/(m + n).
% parameters = [mu, m, n]
parameters = parameters(:).';

if parameters(2) <= 0.0 || parameters(3) <= 0.0
    error('Hill_DN:InvalidParameters', 'Hill-DN requires m > 0 and n > 0.');
end

model.name = 'Hill-DN';
model.parameters = parameters;
model.parameter_names = {'mu', 'm', 'n'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Hill_DN(parameters);
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
n = parameters(3);
out.values = (exp(m .* (lambda - 1.0)) - ...
    exp(n .* (lambda .^ (-1.0) - 1.0))) ./ (m + n);
out.derivatives = (m .* exp(m .* (lambda - 1.0)) + ...
    n .* lambda .^ (-2.0) .* exp(n .* (lambda .^ (-1.0) - 1.0))) ./ (m + n);
end
