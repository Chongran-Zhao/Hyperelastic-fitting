function out = objective(paras, F_UT_list, P_UT_list)
P_UT_fit = get_P_ij_list(1, 1, paras, F_UT_list);

out = P_UT_list - P_UT_fit;
end