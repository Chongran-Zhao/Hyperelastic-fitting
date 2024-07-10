function [paras, lb, ub, num] = array_to_paras(mu, alpha)
num = length(mu);
num_paras = 2*num;
lb = -Inf(num_paras, 1);
ub = Inf(num_paras, 1);
% ub(3) = -1e-16;
% ub(6) = -1e-16;
paras = [mu, alpha];
end
% EOF