function H = get_meas_mat(u_aj_set)
    H = [-1*u_aj_set' ones(size(u_aj_set,2),1)];
end