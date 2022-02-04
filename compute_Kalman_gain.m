function K = compute_Kalman_gain...
    (err_covmat, measurement_mat, measurement_noise_covmat)
    
    P = err_covmat;
    H = measurement_mat;
    R = measurement_noise_covmat;
    
    K = P * H' * inv(H*P*H' + R);
end