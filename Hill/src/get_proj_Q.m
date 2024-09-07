% refer to eqn. 2.8 of Liu, Guan, Zhao & Luo 2024 preprint
% The input C can also be replaced by Gamma, because the transformation
% of E and C is same with Ev and Gamma
function out = get_proj_Q(name, m, n, C)
out = zeros(3,3,3,3);
[V, D] = eig(C);
N1 = V(:,1);
N2 = V(:,2);
N3 = V(:,3);
lambda1 = sqrt(D(1,1));
lambda2 = sqrt(D(2,2));
lambda3 = sqrt(D(3,3));

e1 = get_scale_der(name, m, n, lambda1, 0);
e2 = get_scale_der(name, m, n, lambda2, 0);
e3 = get_scale_der(name, m, n, lambda3, 0);

e1_der = get_scale_der(name, m, n, lambda1, 1);
e2_der = get_scale_der(name, m, n, lambda2, 1);
e3_der = get_scale_der(name, m, n, lambda3, 1);

d1 = e1_der / lambda1;
d2 = e2_der / lambda2;
d3 = e3_der / lambda3;

if lambda1 == lambda2
    theta_12 = d1;
else
    theta_12 = 2.0 * (e1 - e2) / (lambda1*lambda1 - lambda2*lambda2);
end

if lambda2 == lambda3
    theta_23 = d2;
else
    theta_23 = 2.0 * (e2 - e3) / (lambda2*lambda2 - lambda3*lambda3);
end

if lambda1 == lambda3
    theta_13 = d1;
else
    theta_13 = 2.0 * (e1 - e3) / (lambda1*lambda1 - lambda3*lambda3);
end

out = out + d1 .* cross_otimes_1d_to_4d(N1, N1, N1, N1);
out = out + d2 .* cross_otimes_1d_to_4d(N2, N2, N2, N2);
out = out + d3 .* cross_otimes_1d_to_4d(N3, N3, N3, N3);

out = out + theta_12 .* dot_otimes(N1, N2) + theta_12 .* dot_otimes(N2, N1);
out = out + theta_13 .* dot_otimes(N1, N3) + theta_13 .* dot_otimes(N3, N1);
out = out + theta_23 .* dot_otimes(N2, N3) + theta_23 .* dot_otimes(N3, N2);
end