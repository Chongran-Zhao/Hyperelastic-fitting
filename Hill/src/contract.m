% The function is used to calculate contraction between
%    3x3x3x3     fourth-order tensor and 3x3     second-order tensor
% or 3x3x3x3x3x3 sixth-order  tensor and 3x3     second-order tensor
% or 3x3x3x3     fourth-order tensor and 3x3x3x3 fourth-order tensor
function out = contract(A, B)
if length(size(A)) == 4 && length(size(B)) == 2
    out = zeros(3,3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    out(ii,jj) = out(ii,jj) + A(ii,jj,kk,ll) * B(kk,ll);
                end
            end
        end
    end
elseif length(size(B)) == 4 && length(size(A)) == 2
    out = zeros(3,3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    out(kk,ll) = out(kk,ll) + A(ii,jj) * B(ii,jj,kk,ll);
                end
            end
        end
    end
elseif length(size(A)) == 6 && length(size(B)) == 2
    out = zeros(3,3,3,3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    for mm = 1:3
                        for nn = 1:3
                            out(ii,jj,kk,ll) = out(ii,jj,kk,ll) + A(ii,jj,kk,ll,mm,nn) * B(mm,nn);
                        end
                    end
                end
            end
        end
    end
elseif length(size(B)) == 6 && length(size(A)) == 2
    out = zeros(3,3,3,3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    for mm = 1:3
                        for nn = 1:3
                            out(kk,ll,mm,nn) = out(kk,ll,mm,nn) + A(ii,jj) * B(ii,jj,kk,ll,mm,nn);
                        end
                    end
                end
            end
        end
    end
elseif length(size(A)) == 4 && length(size(B)) == 4
    out = zeros(3,3,3,3);
    for ii = 1:3
        for jj = 1:3
            for kk = 1:3
                for ll = 1:3
                    for mm = 1:3
                        for nn = 1:3
                            out(ii,jj,mm,nn) = out(ii,jj,mm,nn) + A(ii,jj,kk,ll) * B(kk,ll,mm,nn);
                        end
                    end
                end
            end
        end
    end
else
    error("ERROR: wrong size of input!");
end
end