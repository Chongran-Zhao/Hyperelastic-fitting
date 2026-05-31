function out = contract(A, B)
if length(size(A)) == 4 && length(size(B)) == 2
    out = zeros(3, 3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    out(ii, jj) = out(ii, jj) + A(ii, jj, kk, ll) .* B(kk, ll);
                end
            end
        end
    end
elseif length(size(B)) == 4 && length(size(A)) == 2
    out = zeros(3, 3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    out(kk, ll) = out(kk, ll) + A(ii, jj) .* B(ii, jj, kk, ll);
                end
            end
        end
    end
elseif length(size(A)) == 2 && length(size(B)) == 2
    out = 0.0;
    for ii = 1:3
        for jj = 1:3
            out = out + A(ii, jj) .* B(ii, jj);
        end
    end
else
    error('contract:UnsupportedContraction', ...
        'Unsupported contraction between arrays of size %s and %s.', ...
        mat2str(size(A)), mat2str(size(B)));
end
end
