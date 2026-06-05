function out = incompressible_constraint(P_iso, F)
%   Apply pressure constraint to incompressible stress.
%
%   out = INCOMPRESSIBLE_CONSTRAINT(P_iso, F) returns the first
%   Piola-Kirchhoff stress
%
%       P = P_iso - p * F^{-T},
%
%   where p is the Lagrange multiplier associated with incompressibility.
%   The pressure is chosen such that the transverse nominal stress vanishes:
%
%       P_33 = 0.
%
%   Therefore,
%
%       p = P_iso(3,3) / F^{-T}(3,3).
%
%   This is appropriate for homogeneous incompressible tests where the
%   third direction is traction-free.

FinvT = F' \ eye(3);
pressure = P_iso(3, 3) ./ FinvT(3, 3);
out = P_iso - pressure .* FinvT;
end