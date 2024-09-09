function out = get_R_square(P_exp, P_pre)
top = 0.0;
bot = 0.0;
for ii = 1:length(P_exp)
top = top + (P_exp(ii) - P_pre(ii))^2;
bot = bot + (P_exp(ii) - mean(P_exp))^2;
end
out = 1.0 - top / bot;
end