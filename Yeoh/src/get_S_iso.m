function out = get_S_iso(xi, C)
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

I1_bar = lambda_bar_1*lambda_bar_1 + lambda_bar_2 * lambda_bar_2 + lambda_bar_3 * lambda_bar_3;
S_bar_1 = Yeoh_dPsi_dlambda(xi, lambda_bar_1, I1_bar) / lambda_bar_1;
S_bar_2 = Yeoh_dPsi_dlambda(xi, lambda_bar_2, I1_bar) / lambda_bar_2;
S_bar_3 = Yeoh_dPsi_dlambda(xi, lambda_bar_3, I1_bar) / lambda_bar_3;
S_bar = tensor_product(S_bar_1, S_bar_2, S_bar_3, N1, N2, N3);
out = detFm0d66 .* DEV(S_bar, C);
end
% EOF