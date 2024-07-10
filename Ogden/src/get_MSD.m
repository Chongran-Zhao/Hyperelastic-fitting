function out = get_MSD(P_exp, P_pre)
out = 0.0;
for ii = 1:length(P_exp)
    out = out + (P_exp(ii)-P_pre(ii)) * (P_exp(ii)-P_pre(ii)) / length(P_exp);
end
end