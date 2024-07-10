% mu and alpha are parameter arrays
% F is the 3x3 deformation gradient at one specific time
% out is a 3x3 2nd PK stress tensor at cooresponding time

function out = get_S_iso(mu, alpha, F)
C = transpose(F)*F;
[V, D] = eig(C);

lambda1 = sqrt(D(1,1));
lambda2 = sqrt(D(2,2));
lambda3 = sqrt(D(3,3));

lambda1_bar = lambda1 / (lambda1*lambda2*lambda3);
lambda2_bar = lambda2 / (lambda1*lambda2*lambda3);
lambda3_bar = lambda3 / (lambda1*lambda2*lambda3);

eig_val_S_iso = zeros(3,1);
eig_val_S_iso(1) = lambda1_bar * Ogden_dPsi_dlambda(mu, alpha, lambda1_bar)...
                -(lambda1_bar * Ogden_dPsi_dlambda(mu, alpha, lambda1_bar)...
                + lambda2_bar * Ogden_dPsi_dlambda(mu, alpha, lambda2_bar)...
                + lambda3_bar * Ogden_dPsi_dlambda(mu, alpha, lambda3_bar)) / 3.0;
eig_val_S_iso(1) = eig_val_S_iso(1) / lambda1^2;

eig_val_S_iso(2) = lambda2_bar * Ogden_dPsi_dlambda(mu, alpha, lambda2_bar)...
                -(lambda1_bar * Ogden_dPsi_dlambda(mu, alpha, lambda1_bar)...
                + lambda2_bar * Ogden_dPsi_dlambda(mu, alpha, lambda2_bar)...
                + lambda3_bar * Ogden_dPsi_dlambda(mu, alpha, lambda3_bar)) / 3.0;
eig_val_S_iso(2) = eig_val_S_iso(2) / lambda2^2;

eig_val_S_iso(3) = lambda3_bar * Ogden_dPsi_dlambda(mu, alpha, lambda3_bar)...
                -(lambda1_bar * Ogden_dPsi_dlambda(mu, alpha, lambda1_bar)...
                + lambda2_bar * Ogden_dPsi_dlambda(mu, alpha, lambda2_bar)...
                + lambda3_bar * Ogden_dPsi_dlambda(mu, alpha, lambda3_bar)) / 3.0;
eig_val_S_iso(3) = eig_val_S_iso(3) / lambda3^2;
out = tensor_product(V, eig_val_S_iso);
end
% EOF