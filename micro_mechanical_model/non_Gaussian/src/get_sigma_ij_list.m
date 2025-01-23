function out = get_sigma_ij_list(ii, jj, paras, F_list)
sigma_list = get_sigma_list(paras, F_list);
out = zeros(size(F_list,3), 1);
out(:) = sigma_list(ii,jj,:);
end
