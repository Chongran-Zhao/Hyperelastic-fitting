function out = Lebedev_quadrature(fun)
srcDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(srcDir));
data = load(fullfile(rootDir, 'data', 'Lebedev.txt'));

phi = data(:, 1);
theta = data(:, 2);
weight = data(:, 3);

out = 0.0;
for ii = 1:length(weight)
    out = out + fun(theta(ii), phi(ii)) .* weight(ii);
end
end
