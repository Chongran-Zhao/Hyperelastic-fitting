function out = get_P_iso(names, paras, F)
C = F' * F;
out = F * get_S_iso(names, paras, C);
end