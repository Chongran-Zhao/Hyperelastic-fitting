function out = get_P_iso(xi, F)
C = F' * F;
out = F * get_S_iso(xi, C);
end