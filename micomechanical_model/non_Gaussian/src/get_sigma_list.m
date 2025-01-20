function out = get_sigma_list(paras, F_list)
out = zeros(size(F_list));
P_list = get_P_list(paras, F_list);
for ii = 1:size(F_list, 3)
    J = det(F_list(:,:,ii));
    out(:,:,ii) = J^(-1.0) .* (P_list(:,:,ii) * F_list(:,:,ii)');
end
end