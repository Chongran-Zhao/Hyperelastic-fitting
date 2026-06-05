function out = dev(input, C)
%   Material deviatoric projection with respect to the C-metric.
%
%   out = input - 1/3 * (input:C) * C^{-1}.
%
%   This corresponds to the material deviatoric operator Dev(.) in
%   Holzapfel (2000), Nonlinear Solid Mechanics, Chapter 6, Eq. (6.90),
%   used in the isochoric second Piola-Kirchhoff stress:
%
%       S_iso = J^(-2/3) * Dev(S_bar).
%
%   Here input is typically the intermediate stress S_bar, and C is the
%   right Cauchy-Green tensor.

out = input - (1.0 / 3.0) .* contract(input, C) .* inv(C);
end