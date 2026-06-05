function out = inv_langevin(lambda_r)
%   Approximate the inverse Langevin function.
%
%   out = INV_LANGEVIN(lambda_r) evaluates the rational approximation
%
%       L^{-1}(x) ≈ x * [3 - 1/5*(6*x^2 + x^4 - 2*x^6)] / (1 - x^2),
%
%   where L(x) = coth(x) - 1/x is the Langevin function.
%
%   This approximation is commonly denoted as L^{-1}_{[3/2]}. It was
%   proposed by Kroger:
%
%       Kroger, M. (2015). Simple, admissible, and accurate approximants of
%       the inverse Langevin and Brillouin functions, relevant for strong
%       polymer deformations and flows. Journal of Non-Newtonian Fluid
%       Mechanics, 223, 77-87.

out = lambda_r ./ (1.0 - lambda_r .* lambda_r) .* ...
    (3.0 - 0.2 .* (6.0 .* lambda_r .* lambda_r + ...
    lambda_r .^ 4.0 - 2.0 .* lambda_r .^ 6.0));
end