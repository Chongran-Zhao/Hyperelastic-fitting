function model = Neo_Hookean(parameters)
% Incompressible Neo-Hookean model.
%
% This model is the Gaussian-network limiting form of rubber elasticity,
% written in terms of the isochoric first invariant
%
%   I1_bar = tr(C_bar) = sum(lambda_bar.^2).
%
% The strain energy density is
%
%   W = 0.5*mu*(I1_bar - 3),
%
% so W(F = I) = 0.
%
% Reference:
%   Treloar, L. R. G. (1943).
%   The elasticity of a network of long-chain molecules-II.
%   Transactions of the Faraday Society, 39, 241-246.
%
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

kin = kinematics(F, 'F');
I1 = sum(kin.lambda_bar .^ 2);
W = 0.5 .* mu .* (I1 - 3.0);
end

function out = S(parameters, F)
mu = parameters(1);

C = F' * F;
kin = kinematics(C, 'C');
S_bar = mu .* eye(3);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
