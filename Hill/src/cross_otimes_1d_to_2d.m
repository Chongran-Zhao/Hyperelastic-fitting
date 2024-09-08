function out = cross_otimes_1d_to_2d(Na, Nb)
out = zeros(3,3);
for ii = 1:3
    for jj = 1:3
                out(ii,jj) = Na(ii) * Nb(jj);
    end
end
end