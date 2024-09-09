function out = get_P_list(names, paras, F_list)
num = size(F_list, 3);
out = zeros(size(F_list));
for ii = 1:num
    P_iso_ii = get_P_iso(names, paras, F_list(:,:,ii));
    out(:,:,ii) = incompressible_constraint(P_iso_ii, F_list(:,:,ii));
end
end