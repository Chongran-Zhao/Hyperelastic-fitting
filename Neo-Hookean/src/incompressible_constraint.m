function out = incompressible_constraint(P_iso, F)
temp = transpose(inv(F));
pressure = P_iso(3,3) / temp(3,3);
out = -pressure .* temp + P_iso;
end