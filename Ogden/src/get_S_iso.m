function out = get_S_iso(num, paras, C)
[V, D] = eig(C);

lambda1 = sqrt(D(1,1));
lambda2 = sqrt(D(2,2));
lambda3 = sqrt(D(3,3));

N1 = V(:,1);
N2 = V(:,2);
N3 = V(:,3);

detF = sqrt(lambda1*lambda2*lambda3);
detFm0d33 = detF^(-1.0/3.0);
detFm0d66 = detF^(-2.0/3.0);

lambda_bar_1 = detFm0d33 * lambda1;
lambda_bar_2 = detFm0d33 * lambda2;
lambda_bar_3 = detFm0d33 * lambda3;

S_bar_1 = 0.0;
S_bar_2 = 0.0;
S_bar_3 = 0.0;
for ii = 1:num
    mu = paras(2*ii-1);
    alpha = paras(2*ii);
    S_bar_1 = S_bar_1 + Ogden_dPsi_dlambda(mu, alpha, lambda_bar_1) / lambda_bar_1;
    S_bar_2 = S_bar_2 + Ogden_dPsi_dlambda(mu, alpha, lambda_bar_2) / lambda_bar_2;
    S_bar_3 = S_bar_3 + Ogden_dPsi_dlambda(mu, alpha, lambda_bar_3) / lambda_bar_3;
end
S_bar = tensor_product(S_bar_1, S_bar_2, S_bar_3, N1, N2, N3);
out = detFm0d66 .* DEV(S_bar, C);
end
% EOF