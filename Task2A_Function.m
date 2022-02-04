function [x, P, lat, long, height, vel] = Task2A_Function(range_data, rangerate_data, x, P, Phi, Q, t, rangeSD, rangerateSD)
    rad = pi/180;
    deg = 180/pi;
    
    %% 2a.d) (Step 3) Use the transition matrix to propagate the state estimates
    x = Phi * x;
    
    %% 2a.e) (Step 4) Then use this to propagate the error covariance matrix
    P = Phi * P * Phi' + Q;
    
    %% 2a.f) Predict the ranges from the approximate user position to each satellite
    % The recursion is resolved by initially computing the range with
    % the Sagnac effect compensation matrix set to the identity matrix, then using this
    % range to compute the Sagnac effect compensation matrix and then recomputing the
    % range.
    [cart_pos, cart_vel] = get_sat_cartesian(range_data, t); % cartesian position and velocity of satellites
    [r_aj, C_Ie] = compute_user2sat_ranges(x(1:3), cart_pos);
    
    %% 2a.g) Compute the line-of-sight unit vector from the approximate user position to each satellite using
    u_aj = calc_line_of_sight_vec(x(1:3), cart_pos, r_aj, C_Ie);
    
    %% 2a.h) Predict the range rates from the approximate user position to each satellite
    rdot_aj = compute_user2sat_range_rates(x(1:3), x(4:6), cart_pos, cart_vel, C_Ie, u_aj);

    %% 2a.i) (Step 5) Compute the measurement matrix
    H = compute_measurement_matrix(u_aj);
 
    %% 2a.j) (Step 6) Compute the measurement noise covariance matrix assuming
    % all pseudorange measurements have an error standard deviation of 10m and all pseudo-range
    % rate measurements have an error standard deviation of 0.05 m/s.
    R = compute_measurement_noise_covmat(range_data, rangeSD, rangerateSD);
    
    %% 2a.k) (Step 7) Compute the Kalman gain matrix
    K = compute_Kalman_gain(P, H, R);
    
    %% 2a.l) (Step 8) Formulate the measurement innovation vector (delta_z)
    dz = compute_measurement_innovation_vec(range_data, rangerate_data, t, x, r_aj, rdot_aj);

    %% 2a.m) (Step 9) Update the state estimates
    x = x + K*dz;
    
    %% 2a.n) (Step 10) Update the error covariance matrix
    P = (eye(8) - K*H)*P;

    %% 2a.o) Convert this Cartesian ECEF position solution to latitude, longitude and height using
    % the Matlab function pv_ECEF_to_NED.m (supplied on Moodle). Note that the
    % Matlab function outputs latitude and longitude in radians. Use the same Matlab
    % function to convert the velocity solution from ECEF resolving axes to north, east and
    % down.
    [lat, long, height, vel] = pv_ECEF_to_NED(x(1:3), x(4:6));
    lat = lat*deg;
    long = long*deg;
end