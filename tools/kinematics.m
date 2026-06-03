function out = kinematics(input, inputType)
% Continuum kinematics for isotropic hyperelastic models.
%
% Polar decomposition:
%   F = R * U
%
% Right Cauchy-Green tensor:
%   C = F' * F = U^2
%
% Spectral decomposition:
%   lambda = [lambda_1, lambda_2, lambda_3]' is the principal stretch
%   vector of U.
%   V = [N_1, N_2, N_3] is the eigenvector matrix from eig(C).
%   M_a = N_a * N_a'
%   C = V * diag(lambda_a^2) * V'
%   C = sum_a lambda_a^2 * M_a
%   U = sum_a lambda_a * M_a
%
% Isochoric split:
%   J = det(F) = prod(lambda_a)
%   U_bar = J^(-1/3) * U
%   C_bar = J^(-2/3) * C
%
% Output convention:
%   U and U_bar are second-order tensors.
%   lambda = [lambda_1, lambda_2, lambda_3]' is the principal stretch
%   vector of U.
%   lambda_bar = [lambda_bar_1, lambda_bar_2, lambda_bar_3]' is the
%   principal stretch vector of U_bar.
%   V stores the principal directions N_a and is not the left stretch tensor.
if nargin < 2
    inputType = 'F';
end

switch upper(char(inputType))
    case 'F'
        F = input;
        C = F' * F;
    case 'C'
        C = input;
    otherwise
        error('kinematics:UnsupportedInputType', ...
            'inputType must be F or C.');
end

[V, D] = eig(C);
lambda = sqrt(max(real(diag(D)), 0.0));
J = prod(lambda);
lambda_bar = J^(-1.0 / 3.0) .* lambda;
U = V * diag(lambda) * V';
U_bar = J^(-1.0 / 3.0) .* U;

out = struct();
out.C = C;
out.V = V;
out.U = U;
out.U_bar = U_bar;
out.J = J;
out.C_bar = J^(-2.0 / 3.0) .* C;
out.lambda = lambda;
out.lambda_bar = lambda_bar;
end
