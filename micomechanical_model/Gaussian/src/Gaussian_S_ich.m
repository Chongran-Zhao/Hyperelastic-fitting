% S_ich = J^(-2/3) * DEV(S_tilde)
function out = Gaussian_S_ich(paras, C)
mu = paras(1);

[lambda_1, lambda_2, lambda_3, N1, N2, N3] = eigen_decomp(C);

J = lambda_1*lambda_2*lambda_3;

tilde_lambda_1 = J^(-1.0/3.0) * lambda_1;
tilde_lambda_2 = J^(-1.0/3.0) * lambda_2;
tilde_lambda_3 = J^(-1.0/3.0) * lambda_3;

S_tilde = 2.0 * mu * ( (tilde_lambda_1+tilde_lambda_2+tilde_lambda_3) ...
        .* tensor_product(1.0/tilde_lambda_1, 1.0/tilde_lambda_2, 1.0/tilde_lambda_3, N1, N2, N3) ...
        + 2.0 .* eye(3) );

out = J^(-2.0/3.0) .* DEV(S_tilde, C);
end