function out = contract(A, B)
%CONTRACT Double contraction between second- and fourth-order tensors.
%
%   This function implements the standard double contraction operation.
%
%   If A and B are both second-order tensors, then
%
%       out = A:B = A_ij B_ij.
%
%   If A is fourth-order and B is second-order, then
%
%       out_ij = A_ijkl B_kl.
%
%   If A is second-order and B is fourth-order, then
%
%       out_kl = A_ij B_ijkl.
%
%   These definitions follow the tensor contraction convention used in
%   Holzapfel (2000), Nonlinear Solid Mechanics, e.g. Eqs. (1.93) and
%   (1.151).

if ismatrix(A) && ismatrix(B)
    validateattributes(A, {'numeric'}, {'size', [3, 3]}, mfilename, 'A');
    validateattributes(B, {'numeric'}, {'size', [3, 3]}, mfilename, 'B');
    out = sum(A .* B, 'all');

elseif ndims(A) == 4 && ismatrix(B)
    validateattributes(A, {'numeric'}, {'size', [3, 3, 3, 3]}, mfilename, 'A');
    validateattributes(B, {'numeric'}, {'size', [3, 3]}, mfilename, 'B');
    out = zeros(3, 3);
    for ii = 1:3
        for jj = 1:3
            out(ii, jj) = sum(squeeze(A(ii, jj, :, :)) .* B, 'all');
        end
    end

elseif ismatrix(A) && ndims(B) == 4
    validateattributes(A, {'numeric'}, {'size', [3, 3]}, mfilename, 'A');
    validateattributes(B, {'numeric'}, {'size', [3, 3, 3, 3]}, mfilename, 'B');
    out = zeros(3, 3);
    for kk = 1:3
        for ll = 1:3
            out(kk, ll) = sum(A .* squeeze(B(:, :, kk, ll)), 'all');
        end
    end

else
    error('contract:UnsupportedContraction', ...
        'Unsupported contraction between arrays of size %s and %s.', ...
        mat2str(size(A)), mat2str(size(B)));
end
end