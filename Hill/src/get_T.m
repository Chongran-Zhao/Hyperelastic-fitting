function out = get_T(name, mu, m, n, C)
out = 2.0 * mu .* get_strain(name, m, n, C);
end