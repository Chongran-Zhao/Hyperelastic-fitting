% Psi = mu * ( tr(U)^2 + 2 * tr(U^2) )
% formula (16) of
% Lin Zhan, Siyu Wang, Shaoxing Qu, Paul Steinmann, Rui Xiao,
% 2023, Journal of the Mechanics and Physics of Solids

function out = Gaussian_energy(mu, C)
e = eig(C);

lambda_1 = sqrt(e(1));
lambda_2 = sqrt(e(2));
lambda_3 = sqrt(e(3));

out = mu * ( (lambda_1 + lambda_2 + lambda_3)^2 + 2.0 * trace(C) );
end
% EOF