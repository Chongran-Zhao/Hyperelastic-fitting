function out = Ogden_dPsi_dlambda(mu, alpha, lambda)
out = 0.0;
for ii = 1:length(mu)
    out = out + mu(ii) * lambda^( alpha(ii)-1.0 );
end
end