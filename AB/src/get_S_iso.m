function out = get_S_iso(paras, C)
% Use devitoric projection
% [V, D] = eig(C);
% lambda1 = sqrt(D(1,1));
% lambda2 = sqrt(D(2,2));
% lambda3 = sqrt(D(3,3));
% 
% detF = sqrt(lambda1*lambda2*lambda3);
% detFm0d33 = detF^(-1.0/3.0);
% detFm0d66 = detF^(-2.0/3.0);
% 
% lambda_bar_1 = detFm0d33 * lambda1;
% lambda_bar_2 = detFm0d33 * lambda2;
% lambda_bar_3 = detFm0d33 * lambda3;
% 
% I_1_bar = lambda_bar_1*lambda_bar_1 + lambda_bar_2*lambda_bar_2 + lambda_bar_3*lambda_bar_3;
% N = paras(2);
% lambda_r_bar = sqrt(I_1_bar / (3.0*N)); 
% mu = paras(1);
% % 2.0 * partial Psi / partial C_bar
% S_bar = sqrt(N/(3.0*I_1_bar)) * mu *inv_Langevin(lambda_r_bar) .* eye(3);
% out = detFm0d66 .* DEV(S_bar, C);

% Use another volumetric energy
[V, D] = eig(C);
lambda1 = sqrt(D(1,1));
lambda2 = sqrt(D(2,2));
lambda3 = sqrt(D(3,3));
I_1 = lambda1*lambda1 + lambda2*lambda2 + lambda3*lambda3;
N = paras(2);
lambda_r = sqrt(I_1 / (3.0*N));
mu = paras(1);
out = sqrt(N/(3.0*I_1)) * mu *inv_Langevin(lambda_r) .* eye(3) - mu * sqrt(N) / 6.0 * inv_Langevin(1.0/sqrt(N)) .* inv(C);
end
% EOF