function out = get_sym_idn_4d(~)
out = zeros(3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(ii,jj,kk,ll) = out(ii,jj,kk,ll) + 0.5 * delta(ii,kk) * delta(jj,ll) + 0.5 * delta(ii,ll) * delta(jj,kk);
            end
        end
    end
end
end