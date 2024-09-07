function out = get_strain(name, m, n, C)
[V, D] = eig(C);
lambda1 = D(1,1);
lambda2 = D(2,2);
lambda3 = D(3,3);
N1 = V(:,1);
N2 = V(:,2);
N3 = V(:,3);
e1 = get_scale_der(name, m, n, lambda1, 0);
e2 = get_scale_der(name, m, n, lambda2, 0);
e3 = get_scale_der(name, m, n, lambda3, 0);
out = tensor_product(e1, e2, e3, N1, N2, N3);
end