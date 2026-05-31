function kinematics = Kinematics(input, inputType)
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
        error('Kinematics:UnsupportedInputType', ...
            'inputType must be F or C.');
end

[V, D] = eig(C);
lambda = sqrt(max(real(diag(D)), 0.0));
J = prod(lambda);

kinematics = struct();
kinematics.C = C;
kinematics.V = V;
kinematics.lambda = lambda;
kinematics.J = J;
kinematics.lambda_bar = J^(-1.0 / 3.0) .* lambda;
kinematics.C_bar = J^(-2.0 / 3.0) .* C;
end
