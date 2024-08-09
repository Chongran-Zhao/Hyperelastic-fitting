function out = cross_otimes_1d_to_4d(Na, Nb, Nc, Nd)
out = zeros(3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(ii,jj,kk,ll) = Na(ii) * Nb(jj) * Nc(kk) * Nd(ll);
            end
        end
    end
end
end