% C is characterized by the right Cauchy-Green deformation tensor
% This function is used to carry out the spectral decomposion

function [lambda_1, lambda_2, lambda_3, N1, N2, N3] = eigen_decomp(C)
[V, D] = eig(C);

lambda_1 = sqrt(D(1,1));
lambda_2 = sqrt(D(2,2));
lambda_3 = sqrt(D(3,3));

N1 = V(:,1);
N2 = V(:,2);
N3 = V(:,3);
end
% EOF