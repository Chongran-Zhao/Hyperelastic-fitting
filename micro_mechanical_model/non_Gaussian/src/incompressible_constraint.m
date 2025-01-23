% The incompressible condition is carried out by
% supposing the P_33 eqauls 0, which is suitable for
% the Uni-bi-axial deformation
function out = incompressible_constraint(P_ich, F)
temp = transpose(inv(F));
pressure = P_ich(3,3) / temp(3,3);
out = -pressure .* temp + P_ich;
end