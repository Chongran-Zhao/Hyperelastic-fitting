function out = kinematics(input, inputType)
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

out = struct();
out.C = C;
out.V = V;
out.lambda = lambda;
out.J = J;
out.lambda_bar = J^(-1.0 / 3.0) .* lambda;
out.C_bar = J^(-2.0 / 3.0) .* C;
end
