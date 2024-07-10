function [mu, alpha] = paras_to_array(paras, num)
mu = paras(1:num);
alpha = paras(num+1:2*num);
end