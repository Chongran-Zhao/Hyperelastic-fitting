function model = Ogden(parameters)
% Incompressible Ogden model with any number of terms.
% W = sum_a mu_a/alpha_a*(sum(lambda_bar.^alpha_a) - 3).
% The number of terms is determined by the number of [mu, alpha] pairs.
% parameters = [mu1, alpha1, mu2, alpha2, ...]
parameters = parameters(:).';

if isempty(parameters) || mod(length(parameters), 2) ~= 0
    error('Ogden:InvalidParameters', ...
        'Ogden parameters must be [mu1, alpha1, mu2, alpha2, ...].');
end

model.name = 'Ogden';
model.parameters = parameters;
model.num_terms = length(parameters) / 2;
model.parameter_names = Ogden_parameter_names(model.num_terms);
model.energy = @(F) Energy(parameters, F);
model.S = @(F) S(parameters, F);
model.P = @(F) P(parameters, F);
model.set_parameters = @(parameters) Ogden(parameters);
end

function names = Ogden_parameter_names(numTerms)
names = cell(1, 2 * numTerms);
for ii = 1:numTerms
    names{2 * ii - 1} = sprintf('mu%d', ii);
    names{2 * ii} = sprintf('alpha%d', ii);
end
end

function W = Energy(parameters, F)
kinematics = Kinematics(F, 'F');
lambda = kinematics.lambda_bar;
numTerms = length(parameters) / 2;
W = 0.0;

for ii = 1:numTerms
    mu = parameters(2 * ii - 1);
    alpha = parameters(2 * ii);
    W = W + mu ./ alpha .* (sum(lambda .^ alpha) - 3.0);
end
end

function out = S(parameters, F)
C = F' * F;
kinematics = Kinematics(C, 'C');
lambda = kinematics.lambda_bar;
numTerms = length(parameters) / 2;
dWdLambda = zeros(3, 1);

for ii = 1:numTerms
    mu = parameters(2 * ii - 1);
    alpha = parameters(2 * ii);
    dWdLambda = dWdLambda + mu .* lambda .^ (alpha - 1.0);
end

S_bar = Spectral_tensor(dWdLambda ./ lambda, kinematics.V);
out = kinematics.J^(-2.0 / 3.0) .* Dev(S_bar, C);
end

function out = P(parameters, F)
out = Incompressible_constraint(F * S(parameters, F), F);
end
