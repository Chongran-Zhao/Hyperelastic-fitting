function out = incompressible_constraint(P_iso, F)
temp = transpose(inv(F));
pressure = P_iso(3, 3) ./ temp(3, 3);
out = P_iso - pressure .* temp;
end
