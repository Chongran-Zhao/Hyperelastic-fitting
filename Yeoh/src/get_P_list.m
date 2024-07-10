function out = get_P_list(xi, Ft)
out = get_P_iso_list(xi, Ft);
% determine the pressure through incompressbility constrain
for ii = 1:size(Ft, 3)
    F = Ft(:,:,ii);
    inv_F_transpose = inv(F');
    p = out(3,3,ii) / inv_F_transpose(3,3);
    out(:,:,ii) = out(:,:,ii) - p .* inv_F_transpose;
end
end
% EOF