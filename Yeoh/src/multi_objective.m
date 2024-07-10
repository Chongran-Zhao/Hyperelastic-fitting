function out = multi_objective(xi, Ft_1, P_exp_1,...
                                      Ft_2, P_exp_2,...
                                      Ft_3, P_exp_3)

P_pre_1 = get_P_ij_list(1, 1, xi, Ft_1);
P_pre_2 = get_P_ij_list(1, 1, xi, Ft_2);
P_pre_3 = get_P_ij_list(1, 1, xi, Ft_3);

out = [];
for ii = 1:length(P_exp_1)
out = [out, P_pre_1(ii) - P_exp_1(ii)];
end
for ii = 1:length(P_exp_2)
out = [out, P_pre_2(ii) - P_exp_2(ii)];
end
for ii = 1:length(P_exp_3)
out = [out, P_pre_3(ii) - P_exp_3(ii)];
end
end