function out = get_P_ij_list(ii, jj, paras, F_list)
P_list = get_P_list(paras, F_list);
out = zeros(size(F_list,3), 1);
out(:) = P_list(ii,jj,:);
end