function out = get_P_list(xi, F_list)
N = size(F_list, 3);
out = zeros(size(F_list));
for ii = 1:N
    P_iso_ii = get_P_iso(xi, F_list(:,:,ii));
    out(:,:,ii) = incompressible_constraint(P_iso_ii, F_list(:,:,ii));
end
end