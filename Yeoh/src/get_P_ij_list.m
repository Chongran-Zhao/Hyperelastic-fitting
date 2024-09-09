function out = get_P_ij_list(ii, jj, xi, F_list)
out = zeros(size(F_list, 3), 1);
P_list = get_P_list(xi, F_list);
out(:) = P_list(ii,jj,:);
end