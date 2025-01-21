% https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
% 1202 point rule, precision 59
function out = Lebedev_quadrature(fun)
out = 0.0;
data = load('Lebedev.txt'); 
phi = data(:,1);
theta = data(:,2);
weight = data(:,3);
for ii = 1:length(weight)
    out = out + fun(theta(ii), phi(ii)) * weight(ii);
end
end