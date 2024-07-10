function out = get_P_ij_list(ii, jj, xi, Ft)
out = zeros(size(Ft, 3), 1);
P_list = get_P_list(xi, Ft);
out(:) = P_list(ii,jj,:);
end