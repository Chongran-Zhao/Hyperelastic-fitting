function out = get_gen_strain(name, m, n, lambda1, lambda2, lambda3, N1, N2, N3)
e1 = get_scale_der(name, m, n, lambda1, 0);
e2 = get_scale_der(name, m, n, lambda2, 0);
e3 = get_scale_der(name, m, n, lambda3, 0);
out = tensor_product(e1, e2, e3, N1, N2, N3);
end