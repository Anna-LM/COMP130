function dz = compute_measurement_innovation_vec(range_data, rangerate_data, t, state_vec, r_aj, rdot_aj)
    range_data = range_data(2:size(range_data,1), :);
    rangerate_data = rangerate_data(2:size(rangerate_data,1), :);
    
    [R1,C1] = find(range_data(:,1) == t);
    [R2,C2] = find(rangerate_data(:,1) == t);    
    
    rho_a = transpose(range_data(R1, 2:size(range_data,2)));
    rhodot_a = transpose(rangerate_data(R2, 2:size(rangerate_data,2)));
    
    r_aj = r_aj';
    rdot_aj = rdot_aj';
    
    clk_offset = state_vec(7) * ones([size(rho_a,1), 1]);
    clk_drift = state_vec(8) * ones([size(rhodot_a,1), 1]);
    
    % disp([rho_a, r_aj, clk_offset]);
    % disp([rhodot_a, rdot_aj, clk_drift]);
    
    dz = [rho_a; rhodot_a] - [r_aj; rdot_aj] - [clk_offset; clk_drift];
end