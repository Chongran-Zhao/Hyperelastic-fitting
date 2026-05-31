function out = Dev(input, C)
out = input - (1.0 / 3.0) .* Contract(input, C) .* inv(C);
end
