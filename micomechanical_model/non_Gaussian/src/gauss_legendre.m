function [points, weights] = gauss_legendre(n)
    % 生成 n 个高斯点和权重
    % 使用 Golub-Welsch 算法
    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    points = diag(D);
    weights = 2 * V(1, :).^2;
end