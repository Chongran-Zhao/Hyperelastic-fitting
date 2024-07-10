function out = get_P_iso_list(mu, alpha, Ft)
out = zeros(size(Ft));
for ii = 1:size(Ft, 3)
    out(:,:,ii) = get_P_iso(mu, alpha, Ft(:,:,ii));
end
end
% EOF