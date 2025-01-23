function val = gauss2D(f, a, b, c, d, n)
% GAUSS2D  在矩形域 [a,b] x [c,d] 上使用 n点 Gauss-Legendre 正交近似
%
% 用法：
%   val = gauss2D(f, a, b, c, d, n)
%
% 输入：
%   f: 函数句柄, 接受标量 (x, y)，返回标量 f(x,y)
%   [a, b]: x 的积分区间
%   [c, d]: y 的积分区间
%   n: 每个方向的 Gauss-Legendre 点数 (越大越精确, 对多项式最高可达 2n-1 次完美积分)
%
% 输出：
%   val: 近似的积分值

    % 1) 一维 Gauss-Legendre 点和权重：在 [-1,1] 上
    [xGL, wGL] = gauss_legendre(n);
    
    % 映射系数
    halfX = (b - a)/2;
    midX  = (b + a)/2;
    
    halfY = (d - c)/2;
    midY  = (d + c)/2;
    
    % 2) 两层求和
    % 先给结果一个初值 0
    val = 0;
    for i = 1:n
        for j = 1:n
            % 映射到 [a,b] x [c,d]
            xx = halfX*xGL(i) + midX;
            yy = halfY*xGL(j) + midY;
            
            % 贡献
            val = val + wGL(i)*wGL(j) * f(xx, yy);
        end
    end
    
    % 不要忘记 dxdy 的映射系数
    val = val * halfX * halfY;
end