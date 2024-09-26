function out = get_P_iso(paras, F)
C = F' * F;
out = F * get_S_iso(paras, C);
end