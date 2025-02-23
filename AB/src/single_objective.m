function out = single_objective(paras, F_list_1, P_exp_1)

P_fit_1 = get_P_ij_list(1, 1, paras, F_list_1);

out = [];
for ii = 1:length(P_fit_1)
out = [out, P_fit_1(ii) - P_exp_1(ii)];
end
end