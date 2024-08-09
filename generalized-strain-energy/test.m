clc;clear;close all
addpath('./src/')
F = rand(3,3);
C = F' * F;
out = get_strain('Curnier-Rakotomanana', 1.0, 1.0, C)
out = get_strain('Darijani-Naghdabadi', 1.0, 1.0, C)