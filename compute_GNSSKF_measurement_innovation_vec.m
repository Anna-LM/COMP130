function dz = compute_GNSSKF_measurement_innovation_vec(pseudo_range, pseudo_rangerate, state_vec, r_aj, rdot_aj)    
    rho_a = transpose(pseudo_range);
    rhodot_a = transpose(pseudo_rangerate);
    
    r_aj = r_aj';
    rdot_aj = rdot_aj';
    
    clk_offset = state_vec(7) * ones([size(rho_a,1), 1]);
    clk_drift = state_vec(8) * ones([size(rhodot_a,1), 1]);
    
    dz = [rho_a; rhodot_a] - [r_aj; rdot_aj] - [clk_offset; clk_drift];
end