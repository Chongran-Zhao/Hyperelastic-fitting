function out = DEV(C)
out = get_sym_idn_4d();
out = out - cross_otimes_2d_to_4d(inv(C), C) ./ 3.0;
end