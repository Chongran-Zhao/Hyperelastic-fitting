function out = get_P_ij_list(ii, jj, mu, alpha, Ft)
out = zeros(size(Ft, 3), 1);
P_list = get_P_list(mu, alpha, Ft);
out(:) = P_list(ii,jj,:);
end