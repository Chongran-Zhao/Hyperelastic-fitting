function out = inv_Langevin(lambda_r)
out = lambda_r / (1.0 - lambda_r*lambda_r) * (3.0 - 0.2*(6.0*lambda_r*lambda_r + lambda_r^4.0 - 2.0*lambda_r^6.0));
end