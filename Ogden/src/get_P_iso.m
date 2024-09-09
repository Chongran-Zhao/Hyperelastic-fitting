function out = get_P_iso(num, paras, F)
C = F' * F;
out = F * get_S_iso(num, paras, C);
end