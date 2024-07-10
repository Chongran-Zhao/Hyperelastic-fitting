% D is a 1x3 vector containing three eigenvalues
% V is a 3x3 matrix containing three cooresponidng eigenvectors
function out = tensor_product(V, D)
out = zeros(3,3);
out = out + D(1) .* kron(V(:,1), V(:,1)');
out = out + D(2) .* kron(V(:,2), V(:,2)');
out = out + D(3) .* kron(V(:,3), V(:,3)');
end