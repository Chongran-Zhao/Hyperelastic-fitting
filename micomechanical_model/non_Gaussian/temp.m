clc;clear;close all

f = @(theta, phi)  cos(theta)*cos(theta);
g = @(theta, phi)  sin(theta)*sin(theta)*cos(phi)*cos(phi);
h = @(theta, phi)  sin(theta)*sin(theta)*sin(phi)*sin(phi);

out = Lebedev_quadrature(f)
out = Lebedev_quadrature(g)
out = Lebedev_quadrature(h)