% P_ich = F * S_ich
% C = F^T * F
function out = get_P_ich(paras, F)
out = F * non_Gaussian_S_ich(paras, F' * F);
end