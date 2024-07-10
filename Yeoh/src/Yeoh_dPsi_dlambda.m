function out = Yeoh_dPsi_dlambda(xi, lambda, I1)
C1 = xi(1);
C2 = xi(2);
C3 = xi(3);
out = 2.0 * lambda * (C1 + 2.0*C2*(I1-3.0) + 3*C3*(I1-3.0)*(I1-3.0));
end