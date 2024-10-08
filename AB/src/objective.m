function out = objective(paras, F_list_1, P_exp_1,...
                                F_list_2, P_exp_2,...
                                F_list_3, P_exp_3)

P_fit_1 = get_P_ij_list(1, 1, paras, F_list_1);
P_fit_2 = get_P_ij_list(1, 1, paras, F_list_2);
P_fit_3 = get_P_ij_list(1, 1, paras, F_list_3);

out = [];
for ii = 1:length(P_fit_1)
out = [out, P_fit_1(ii) - P_exp_1(ii)];
end
for ii = 1:length(P_exp_2)
out = [out, P_fit_2(ii) - P_exp_2(ii)];
end
for ii = 1:length(P_exp_3)
out = [out, P_fit_3(ii) - P_exp_3(ii)];
end
end