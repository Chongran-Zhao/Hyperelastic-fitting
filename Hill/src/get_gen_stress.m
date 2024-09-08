function out = get_gen_stress(mu, gen_strain)
    out = 2.0 * mu .* gen_strain;
end