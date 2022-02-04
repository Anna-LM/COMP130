function dz = get_dz_vec(pseud_ranges, range, clock_offset)
    rho = pseud_ranges;
    r_a = range;

    total_offset = clock_offset*ones(size(rho,2),1);

    dz = transpose(rho) - transpose(r_a) - total_offset;
    
end