function dz = get_meas_innovation_vec(pseudo_ranges, sat2user, offset)
    rho = pseudo_ranges;
    r = sat2user;
    
    dz = rho' - r' - offset*ones(size(rho,2),1);
end