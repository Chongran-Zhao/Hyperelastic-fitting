function out = lebedev_quadrature(fun)
%   Lebedev quadrature on the unit sphere.
%
%   out = LEBEDEV_QUADRATURE(fun) evaluates the spherical average
%
%       out = sum_i weight_i * fun(theta_i, phi_i),
%
%   where the weights are normalized such that sum_i weight_i = 1. Thus the
%   returned value approximates
%
%       (1/S) * integral_{S^2} fun(theta, phi) ds.
%
%   The data file is the 110-point, precision-17 Lebedev rule from
%
%       https://people.math.sc.edu/Burkardt/datasets/sphere_lebedev_rule/lebedev_017.txt
%
%   The first column is the azimuthal angle phi in degrees, the second
%   column is the polar angle theta in degrees, and the third column is the
%   normalized quadrature weight.
%
%   The function handle fun must accept inputs as fun(theta, phi).

srcDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(srcDir);
data = load(fullfile(rootDir, 'data', 'Lebedev.txt'));

phi = deg2rad(data(:, 1));
theta = deg2rad(data(:, 2));
weight = data(:, 3);

out = fun(theta(1), phi(1)) .* weight(1);
for ii = 2:length(weight)
    out = out + fun(theta(ii), phi(ii)) .* weight(ii);
end
end