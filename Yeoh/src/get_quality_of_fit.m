function out = get_quality_of_fit(P_pre, P_exp)
out = 0.0;
for ii = 1:length(P_exp)
out = out + (P_pre(ii) - P_exp(ii))^2 / P_exp(ii);
end