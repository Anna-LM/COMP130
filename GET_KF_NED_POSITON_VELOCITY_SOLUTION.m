function results = GET_KF_NED_POSITON_VELOCITY_SOLUTION(...
    DR_calculated_data, GNSS_KF_data, tau, PSD)
    
    Define_Constants;

    time = DR_calculated_data(:,1);
    
    DR_L = DR_calculated_data(:,2)*rad;
    DR_lamda = DR_calculated_data(:,3)*rad;
    DR_vN = DR_calculated_data(:,4);
    DR_vE = DR_calculated_data(:,5);
    
    GNSS_L = GNSS_KF_data(:,2)*rad;
    GNSS_lamda = GNSS_KF_data(:,3)*rad;
    GNSS_h = GNSS_KF_data(:,4);
    GNSS_vN = GNSS_KF_data(:,5);
    GNSS_vE = GNSS_KF_data(:,6);
    
    results = zeros(size(time,1), 4);
    
    SDGr = 5.2;  % signal in space, ionosphere, troposphere, multipath error
    SDGv = 0.05; % wheel speed sensor noise
    
    %% measurement matrix
    H = zeros(4);
    H(1,3) = -1; H(2,4) = -1; H(3,1) = -1; H(4,2) = -1;
    
    %% Initialise error state vector
    x = zeros(4,1);
    
    %% Initialise error covariance matrix
    [RN, RE] = Radii_of_curvature(DR_L(1));
    P = [1; 1; 10^2/(RN+GNSS_h(1))^2; 10^2/((RE+GNSS_h(1))^2*cos(DR_L(1)))];
    P = diag(P);
    
    for i = 1:size(time,1)
        Psi = compute_transition_mat(tau, DR_L(i), GNSS_h(i));
        Q = compute_sys_noise_covmat(tau, PSD, DR_L(i), GNSS_h(i));
        x = Psi * x;
        P = Psi * P * Psi' + Q;
        R = integrated_measurement_noise_covmat(SDGr, SDGv, DR_L(i), GNSS_h(i));
        K = P * H' * (H*P*H' + R);
        z = [GNSS_L(i) - DR_L(i);
             GNSS_lamda(i) - DR_lamda(i);
             GNSS_vN(i) - DR_vN(i);
             GNSS_vE(i) - DR_vE(i)];
        dz = z - H*x;
        x = x + K*dz;
        
        results(i,1) = (DR_L(i) - x(3))*deg;
        results(i,2) = (DR_lamda(i) - x(4))*deg;
        results(i,3) = DR_vN(i) - x(1);
        results(i,4) = DR_vE(i) - x(2);
        
        P = (eye(4) - K*H) * P;
    end
    
    results = [time results];
end