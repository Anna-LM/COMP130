function H = compute_GNSSKF_measurement_matrix(LoS_vec)
    u_aj = -1*LoS_vec';
    Z_mat = zeros(size(u_aj));          % zero matrix with same dimensions as u_aj set
    Z_vec = zeros([size(u_aj,1), 1]);   % zero vector with same no. of rows as u_aj set
    I_vec = ones([size(u_aj,1), 1]);    % unit vector with same no. of rows as u_aj set
    
    H = [u_aj Z_mat I_vec Z_vec;
         Z_mat u_aj Z_vec I_vec];
end