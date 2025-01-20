function out = inv_Langevin(x)
if isscalar(x)
    out = x / (1.0 - x^2) * (3.0 - 0.2 * (6.0 * x^2 + x^4 - 2.0 * x^6));
else
    out = x ./ (ones(size(x)) - x.^2) .* (3.0 * ones(size(x)) - 0.2 * (6.0 * x.^2 + x.^4 - 2.0 * x.^6));
end
end