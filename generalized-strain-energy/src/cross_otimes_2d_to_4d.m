function out = cross_otimes_2d_to_4d(M1, M2)
out = zeros(3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(ii,jj,kk,ll) = M1(ii,jj) * M2(kk,ll);
            end
        end
    end
end
end