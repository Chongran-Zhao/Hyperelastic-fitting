function out = non_Gaussian_S_ich(paras, C)
    mu = paras(1);
    N = paras(2);

    [lambda_1, lambda_2, lambda_3, N1, N2, N3] = eigen_decomp(C);

    J = lambda_1 * lambda_2 * lambda_3;

    tilde_lambda_1 = J^(-1.0/3.0) * lambda_1;
    tilde_lambda_2 = J^(-1.0/3.0) * lambda_2;
    tilde_lambda_3 = J^(-1.0/3.0) * lambda_3;

    lambda_i = @(theta, phi) tilde_lambda_1 .* cos(theta).^2 + ...
        tilde_lambda_2 .* sin(theta).^2 .* cos(phi).^2 + ...
        tilde_lambda_3 .* sin(theta).^2 .* sin(phi).^2;

    beta = @(theta, phi) inv_Langevin(lambda_i(theta, phi) ./ sqrt(N));

    temp_1 = @(theta, phi) beta(theta, phi) .* cos(theta) .* cos(theta) .* sin(theta);
    temp_2 = @(theta, phi) beta(theta, phi) .* sin(theta) .* sin(theta) .* cos(phi) .* cos(phi) .* sin(theta);
    temp_3 = @(theta, phi) beta(theta, phi) .* sin(theta) .* sin(theta) .* sin(phi) .* sin(phi) .* sin(theta);

    % 使用高斯积分计算
    S_tilde_1 = mu * sqrt(N) / tilde_lambda_1 * gauss2D(temp_1, 0, pi, 0, 2*pi, 10);
    S_tilde_2 = mu * sqrt(N) / tilde_lambda_2 * gauss2D(temp_2, 0, pi, 0, 2*pi, 10);
    S_tilde_3 = mu * sqrt(N) / tilde_lambda_3 * gauss2D(temp_3, 0, pi, 0, 2*pi, 10);

    S_tilde = tensor_product(S_tilde_1, S_tilde_2, S_tilde_3, N1, N2, N3);
    out = J^(-2.0/3.0) .* DEV(S_tilde, C);
end