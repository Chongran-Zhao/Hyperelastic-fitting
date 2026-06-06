function model = Zhan_Gaussian(parameters)
% Zhan Gaussian chain-network model.
%
% The model follows the micro-macro transition proposed by Zhan et al.
% using the isochoric right stretch tensor U_bar. In the principal basis,
%
%   U_bar = sum_a lambda_bar_a * M_a.
%
% For Gaussian-chain statistics, the directional chain stretch is averaged
% analytically over the unit sphere, giving the closed-form energy
%
%   W = mu*((tr(U_bar))^2 + 2*tr(U_bar^2) - 15).
%
% The constant -15 makes W(F = I) = 0 because U_bar = I in the reference
% configuration.
%
% Reference:
%   Zhan, L., Wang, S., Qu, S., Steinmann, P., & Xiao, R. (2023).
%   A new micro-macro transition for hyperelastic materials.
%   Journal of the Mechanics and Physics of Solids, 171, 105156.
%
% parameters = [mu]
parameters = parameters(:).';

model.name = 'Zhan-Gaussian';
model.parameters = parameters;
model.parameter_names = {'mu'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Zhan_Gaussian(parameters);
end

function W = Energy(parameters, F)
kin = kinematics(F, 'F');
mu = parameters(1);
lambda_bar = kin.lambda_bar;
W = mu .* (sum(lambda_bar) .^ 2 + 2.0 .* sum(lambda_bar .^ 2) - 15.0);
end

function out = S(parameters, F)
C = F' * F;
kin = kinematics(C, 'C');
mu = parameters(1);
lambda_bar = kin.lambda_bar;

values = 2.0 .* mu .* ((sum(lambda_bar) ./ lambda_bar) + 2.0);
S_bar = spectral_tensor(values, kin.V);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
