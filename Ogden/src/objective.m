function out = objective(paras, Ft, P_exp, num)
[mu, alpha] = paras_to_array(paras, num);

P_pre = get_P_ij_list(1, 1, mu, alpha, Ft);
% out = [];
% for ii = 1:length(P_exp)
% out = [out, P_exp(ii) - P_pre(ii)];
% end
out = 0.0;
for ii = 1:length(P_exp)
out = out + (P_exp(ii) - P_pre(ii))^2 / length(P_exp);
end
end