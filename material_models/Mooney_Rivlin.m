function model = Mooney_Rivlin(parameters)
% Incompressible Mooney-Rivlin model.
%
% The model is written in terms of the isochoric invariants
%
%   I1_bar = tr(C_bar),
%   I2_bar = 1/2*(I1_bar^2 - tr(C_bar^2)),
%
% which are evaluated from the isochoric principal stretches lambda_bar.
% The strain energy density is
%
%   W = C1*(I1_bar - 3) + C2*(I2_bar - 3),
%
% so W(F = I) = 0.
%
% References:
%   Mooney, M. (1940).
%   A theory of large elastic deformation.
%   Journal of Applied Physics, 11(9), 582-592.
%
%   Rivlin, R. S. (1948).
%   Large elastic deformations of isotropic materials IV. Further
%   developments of the general theory.
%   Philosophical Transactions of the Royal Society of London. Series A,
%   Mathematical and Physical Sciences, 241(835), 379-397.
%
% parameters = [C1, C2]
parameters = parameters(:).';

model.name = 'Mooney-Rivlin';
model.parameters = parameters;
model.parameter_names = {'C1', 'C2'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Mooney_Rivlin(parameters);
end

function W = Energy(parameters, F)
C1 = parameters(1);
C2 = parameters(2);

kin = kinematics(F, 'F');
lambda = kin.lambda_bar;
I1 = sum(lambda .^ 2);
I2 = lambda(1)^2 * lambda(2)^2 + ...
    lambda(2)^2 * lambda(3)^2 + ...
    lambda(1)^2 * lambda(3)^2;

W = C1 .* (I1 - 3.0) + C2 .* (I2 - 3.0);
end

function out = S(parameters, F)
C1 = parameters(1);
C2 = parameters(2);

C = F' * F;
kin = kinematics(C, 'C');
lambda = kin.lambda_bar;

values = zeros(3, 1);
for ii = 1:3
    other = setdiff(1:3, ii);
    values(ii) = 2.0 .* (C1 + C2 .* sum(lambda(other) .^ 2));
end

S_bar = spectral_tensor(values, kin.V);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
