function out = get_P_list(paras, F_list)
out = zeros(size(F_list));
for ii = 1:size(F_list, 3)
    P_ich_ii = get_P_ich(paras, F_list(:,:,ii));
    out(:,:,ii) = incompressible_constraint(P_ich_ii, F_list(:,:,ii));
end
end