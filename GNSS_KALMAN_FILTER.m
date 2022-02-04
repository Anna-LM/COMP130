function [GNSSKF_cartpos, GNSSKF_cartvel] = GNSS_KALMAN_FILTER(range_data, rangerate_data, x_init, P_init, Phi, Q, rangeSD, rangerateSD, sat2user_range, sat2user_rangerate, LoS)
    t = range_data(2:size(range_data,1),1);
    sat = range_data(1,2:size(range_data,2));
    range_data = range_data(2:size(range_data,1),2:size(range_data,2));
    rangerate_data = rangerate_data(2:size(rangerate_data,1),2:size(rangerate_data,2));
    x = x_init;
    P = P_init;
    
    GNSSKF_cartpos = zeros(size(t,1),3);
    GNSSKF_cartvel = GNSSKF_cartpos;
    
    for i = 1:size(t,1)
        %% 2a.d) (Step 3) Use the transition matrix to propagate the state estimates
        x = Phi * x;

        %% 2a.e) (Step 4) Then use this to propagate the error covariance matrix
        P = Phi * P * Phi' + Q;

        %% 2a.f) Predict the ranges from the approximate user position to each satellite
        r_aj = sat2user_range(i,:);

        %% 2a.g) Compute the line-of-sight unit vector from the approximate user position to each satellite
        u_aj = LoS(:,:,i);

        %% 2a.h) Predict the range rates from the approximate user position to each satellite
        rdot_aj = sat2user_rangerate(i,:);

        %% 2a.i) (Step 5) Compute the measurement matrix
        H = compute_GNSSKF_measurement_matrix(u_aj);

        %% 2a.j) (Step 6) Compute the measurement noise covariance matrix
        R = compute_measurement_noise_covmat(sat, rangeSD, rangerateSD);

        %% 2a.k) (Step 7) Compute the Kalman gain matrix
        K = compute_Kalman_gain(P, H, R);

        %% 2a.l) (Step 8) Formulate the measurement innovation vector (delta_z)
        dz = compute_GNSSKF_measurement_innovation_vec(range_data(i,:), rangerate_data(i,:), x, r_aj, rdot_aj);

        %% 2a.m) (Step 9) Update the state estimates
        x = x + K*dz;

        %% 2a.n) (Step 10) Update the error covariance matrix
        P = (eye(8) - K*H)*P;
        
        %% store results
        GNSSKF_cartpos(i,:) = transpose(x(1:3));
        GNSSKF_cartvel(i,:) = transpose(x(4:6));
    end
end