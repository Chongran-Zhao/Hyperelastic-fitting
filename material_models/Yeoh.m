function model = Yeoh(parameters)
% Yeoh model driven by the first isochoric invariant.
% I1_bar = sum(lambda_bar.^2).
% W = C1*(I1_bar - 3) + C2*(I1_bar - 3)^2 + C3*(I1_bar - 3)^3.
% parameters = [C1, C2, C3]
parameters = parameters(:).';

model.name = 'Yeoh';
model.parameters = parameters;
model.parameter_names = {'C1', 'C2', 'C3'};
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Yeoh(parameters);
end

function W = Energy(parameters, F)
C1 = parameters(1);
C2 = parameters(2);
C3 = parameters(3);

kin = kinematics(F, 'F');
I1 = sum(kin.lambda_bar .^ 2);

W = C1 .* (I1 - 3.0) + ...
    C2 .* (I1 - 3.0) .^ 2 + ...
    C3 .* (I1 - 3.0) .^ 3;
end

function out = S(parameters, F)
C1 = parameters(1);
C2 = parameters(2);
C3 = parameters(3);

C = F' * F;
kin = kinematics(C, 'C');
lambda = kin.lambda_bar;
I1 = sum(lambda .^ 2);

dWdLambda = 2.0 .* lambda .* ...
    (C1 + 2.0 .* C2 .* (I1 - 3.0) + ...
    3.0 .* C3 .* (I1 - 3.0) .^ 2);

S_bar = spectral_tensor(dWdLambda ./ lambda, kin.V);
out = kin.J^(-2.0 / 3.0) .* dev(S_bar, C);
end

function out = P(parameters, F)
out = incompressible_constraint(F * S(parameters, F), F);
end
