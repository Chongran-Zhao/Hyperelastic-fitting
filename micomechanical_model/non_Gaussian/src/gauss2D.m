function I = gauss2D(f, a, b, c, d, n)
    % f: 被积函数
    % a, b: x 的积分区间
    % c, d: y 的积分区间
    % n: 高斯点的数量

    % 生成高斯点和权重
    [points_x, weights_x] = gauss_legendre(n);
    [points_y, weights_y] = gauss_legendre(n);

    % 映射到实际积分区间
    points_x = 0.5 * (points_x + 1) * (b - a) + a;
    points_y = 0.5 * (points_y + 1) * (d - c) + c;

    % 计算加权和
    I = 0;
    for i = 1:n
        for j = 1:n
            x = points_x(i);
            y = points_y(j);
            weight = 0.5 * (b - a) * weights_x(i) * 0.5 * (d - c) * weights_y(j);
            I = I + weight * f(x, y);
        end
    end
end