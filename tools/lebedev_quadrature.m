function out = lebedev_quadrature(fun)
% Lebedev quadrature on the unit sphere.
%
% The data file is the 110-point, precision-17 Lebedev rule from:
%   https://people.math.sc.edu/Burkardt/datasets/sphere_lebedev_rule/lebedev_017.txt
% Its first two columns are angles in degrees and the weights are
% normalized so that sum(weight) = 1.
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
