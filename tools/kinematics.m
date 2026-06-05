function out = kinematics(input, inputType)
%   Continuum kinematics for isotropic hyperelastic models.
%
%   out = KINEMATICS(input, inputType) computes the basic kinematic
%   quantities used in isotropic hyperelastic constitutive models. The input
%   can be either the deformation gradient F or the right Cauchy-Green tensor
%   C, specified by inputType = 'F' or inputType = 'C'.
%
%   Polar decomposition:
%
%       F = R * U.
%
%   Right Cauchy-Green tensor:
%
%       C = F' * F = U^2.
%
%   Spectral decomposition:
%
%       C = sum_a lambda_a^2 * M_a,
%       U = sum_a lambda_a   * M_a,
%       M_a = N_a * N_a'.
%
%   Here lambda = [lambda_1, lambda_2, lambda_3]' is the principal stretch
%   vector, and V = [N_1, N_2, N_3] stores the corresponding principal
%   directions returned by eig(C). The ordering of the principal directions
%   follows MATLAB's eig output.
%
%   Isochoric split:
%
%       J = det(F) = sqrt(det(C)) = prod(lambda_a),
%       U_bar = J^(-1/3) * U,
%       C_bar = J^(-2/3) * C,
%       lambda_bar_a = J^(-1/3) * lambda_a.
%
%   Output fields:
%
%       out.C          right Cauchy-Green tensor
%       out.V          principal direction matrix [N_1, N_2, N_3]
%       out.U          right stretch tensor
%       out.U_bar      isochoric right stretch tensor
%       out.C_bar      isochoric right Cauchy-Green tensor
%       out.J          volume ratio
%       out.lambda     principal stretches of U
%       out.lambda_bar isochoric principal stretches of U_bar
%
%   Note that V denotes the eigenvector matrix of C and is not the left
%   stretch tensor.

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
