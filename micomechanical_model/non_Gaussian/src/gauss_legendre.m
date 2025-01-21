function [x, w] = gauss_legendre(n)
% 生成 n点 Gauss-Legendre 正交的节点 x, 权重 w, 适用于 [-1,1]
%
% 简化实现，仅供示例参考，可在 File Exchange 找更健壮版本

    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));  % Golub-Welsch
    T = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(T);
    x = diag(D);   % 节点
    w = 2*(V(1,:).^2); % 权重
    
    % 排序 (升序)
    [x, idx] = sort(x);
    w = w(idx);
end