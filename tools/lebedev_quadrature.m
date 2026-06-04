function out = lebedev_quadrature(fun)
srcDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(srcDir);
data = load(fullfile(rootDir, 'data', 'Lebedev.txt'));

phi = data(:, 1);
theta = data(:, 2);
weight = data(:, 3);

out = fun(theta(1), phi(1)) .* weight(1);
for ii = 2:length(weight)
    out = out + fun(theta(ii), phi(ii)) .* weight(ii);
end
end
