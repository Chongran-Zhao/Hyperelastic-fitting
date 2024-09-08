% refer to eqn. 2.10 of Liu, Guan, Zhao & Luo 2024 preprint
function out = dot_otimes(Na, Nb)
out = 0.5 .* cross_otimes_1d_to_4d(Na, Nb, Na, Nb) + 0.5 .* cross_otimes_1d_to_4d(Na, Nb, Nb, Na);
end