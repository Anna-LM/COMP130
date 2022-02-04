function results = GET_INTEGRATED_KF_SOLUTION(...
    DR_calculated_data, gyro_heading, DR_original_data, GNSS_Kalman_Data, tau_s, wheel_PSD,...
    gyro_meas_noise, gyro_noise, gyro_bias, mag_noiseSD)
    
    Define_Constants;

    time = DR_calculated_data(:,1);
    
    DR_L = DR_calculated_data(:,2)*rad;
    DR_lamda = DR_calculated_data(:,3)*rad;
    DR_vN = DR_calculated_data(:,4);
    DR_vE = DR_calculated_data(:,5);
    
    GNSS_L = GNSS_Kalman_Data(:,2)*rad;
    GNSS_lamda = GNSS_Kalman_Data(:,3)*rad;
    GNSS_h = GNSS_Kalman_Data(:,4);
    GNSS_vN = GNSS_Kalman_Data(:,5);
    GNSS_vE = GNSS_Kalman_Data(:,6);
    
    gyro = gyro_heading;
    compass = DR_original_data(:,7)*rad;
    
    results = zeros(size(time,1), 5);
    
    SDGr = 5.2;  % signal in space, ionosphere, troposphere, multipath error
    SDGv = 0.05; % wheel speed sensor noise
    
    %% measurement matrix
    H = zeros(6);
    H(1,3) = -1; H(2,4) = -1; H(3,1) = -1; H(4,2) = -1; H(5,5) = -1;
    
    %% Initialise error state vector
    x = zeros(6,1);
    
    %% Initialise error covariance matrix
    [RN, RE] = Radii_of_curvature(DR_L(1));
    P = [0.01; 0.01; 10^2/(RN+GNSS_h(1))^2; 10^2/((RE+GNSS_h(1))^2*cos(DR_L(1))); 3e-6; rad^2];
    P = diag(P);
    
    for i = 1:size(time,1)
        Psi = compute_transition_matrix(tau_s, DR_L(i), GNSS_h(i));
        Q = compute_sys_noise_covmatrix(tau_s, wheel_PSD, DR_L(i), GNSS_h(i), gyro_noise, gyro_bias);
        x = Psi * x;
        P = Psi * P * Psi' + Q;
        R = integrated_measurement_noise_covmatrix(SDGr, SDGv, DR_L(i), GNSS_h(i), mag_noiseSD);
        K = P * H' * (H*P*H' + R);
        z = [GNSS_L(i) - DR_L(i);
             GNSS_lamda(i) - DR_lamda(i);
             GNSS_vN(i) - DR_vN(i);
             GNSS_vE(i) - DR_vE(i);
             compass(i) - gyro(i);
             0];
        dz = z - H*x;
        x = x + K*dz;
        
        results(i,1) = (DR_L(i) - x(3))*deg;
        results(i,2) = (DR_lamda(i) - x(4))*deg;
        results(i,3) = DR_vN(i) - x(1);
        results(i,4) = DR_vE(i) - x(2);
        results(i,5) = (gyro(i) - x(5))*deg;
        
        P = (eye(6) - K*H) * P;
    end
    
    results = [time results];
end