function out = spectral_tensor(values, V)
out = zeros(3, 3);
for ii = 1:3
    Ni = V(:, ii);
    out = out + values(ii) .* (Ni * Ni');
end
end
