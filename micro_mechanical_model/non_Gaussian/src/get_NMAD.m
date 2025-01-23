function out = get_NMAD(P_pre, P_exp)
aa = 0.0;
bb = 0.0;
cc = 0.0;
for ii = 1:length(P_pre)
    aa = aa + abs(P_pre(ii) - P_exp(ii));
    bb = bb + abs(P_exp(ii));
    cc = cc + abs(P_pre(ii));
end
aa = aa / length(P_pre);
bb = bb / length(P_exp);
cc = cc / length(P_pre);
out = 100 * aa / max(bb, cc);
end