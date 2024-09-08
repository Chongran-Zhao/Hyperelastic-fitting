function out = get_P_list(names, paras, F_list)
num = size(F_list, 3);
out = zeros(size(F_list));
for ii = 1:num
    F = F_list(:,:,ii);
    P_iso_ii = get_P_iso(names, paras, F);
    out(:,:,ii) = incompressible_constraint(P_iso_ii, F);
end
end