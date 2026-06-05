function out = spectral_tensor(values, V)
%   Reconstruct a second-order tensor from spectral data.
%
%   out = SPECTRAL_TENSOR(values, V) returns
%
%       out = sum_a values(a) * N_a * N_a',
%
%   where values(a) are the principal values and V = [N_1, N_2, N_3]
%   contains the corresponding principal directions as its columns.
%
%   This is used for isotropic tensor functions written in spectral form,
%   for example
%
%       A = sum_a A_a M_a,    M_a = N_a * N_a'.
%
%   The ordering of values must be consistent with the ordering of the
%   eigenvectors stored in V.

out = zeros(3, 3);
for ii = 1:3
    Ni = V(:, ii);
    out = out + values(ii) .* (Ni * Ni');
end
end