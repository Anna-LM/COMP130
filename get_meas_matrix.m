function H_g = get_meas_matrix(u_aj)
    fourth_col = ones(size(u_aj,2),1);
    H_g = [-1*transpose(u_aj) fourth_col];
end