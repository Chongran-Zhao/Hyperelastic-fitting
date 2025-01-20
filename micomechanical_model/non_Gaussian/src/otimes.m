% N1 and N2 are abstract column vector with the size of 3x1
function out = otimes(N1, N2)
out = kron(N2, N2');
end