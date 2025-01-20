function out = non_Gaussian_energy(paras, C)
mu = paras(1);
N = paras(2);
e = eig(C);

lambda_1 = sqrt(e(1));
lambda_2 = sqrt(e(2));
lambda_3 = sqrt(e(3));

lambda_i = @(theta, phi) lambda_1 .* cos(theta).^2 + ...
    lambda_2 .* sin(theta).^2 .* cos(phi).^2 + ...
    lambda_3 .* sin(theta).^2 .* sin(phi).^2;

beta = @(theta, phi) inv_Langevin(lambda_i(theta, phi) ./ sqrt(N));

Psi = @(theta, phi) mu * N .* (lambda_i(theta, phi) .* beta(theta, phi) ./ sqrt(N) + ...
    log(beta(theta, phi) ./ sinh(beta(theta, phi)))) .* sin(theta);

out = integral2(Psi, 0, pi, 0, 2*pi, 'AbsTol', 1e-12, 'RelTol', 1e-9);
end