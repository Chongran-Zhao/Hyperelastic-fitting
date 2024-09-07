function out = get_S_iso(name, mu, m, n, C)
C_bar = C ./ (sqrt(det(C)))^(2.0/3.0);
S_bar = contract(get_T(name, mu, m, n, C_bar), get_proj_Q(name, m, n, C_bar));
P = get_proj_P(C);
out = contract(P, S_bar) ./ (sqrt(det(C)))^(2.0/3.0);
end