function out = DEV(input, C)
out = input - (1.0/3.0) * contract(input, C) .* inv(C);
end