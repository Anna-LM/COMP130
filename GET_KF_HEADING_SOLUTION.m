function heading = GET_KF_HEADING_SOLUTION(...
    DR_data, tau, gyro_meas_noise, mag_noise_SD, gyro_noise, gyro_bias, gyro_heading)
    
    Define_Constants;
    
    %% extract data from dead reckoning data
    time = DR_data(:,1);
    gyro = gyro_heading;
    compass = DR_data(:,7)*rad;
    
    heading = zeros(size(time));
    
    %% Transition matrix
    Psi = [1 tau;
           0 1];
       
    %% System noise covariance matrix
    Q = [gyro_noise*tau+gyro_bias/3*tau^3 0.5*gyro_bias*tau^2;
         0.5*gyro_bias*tau^2 gyro_bias*tau];
    
    %% Initialise error state vector
    x = [0;     % gyro error
         0];    % gyro bias
     
    %% Initialise error covariance matrix
    %P = [1; rad^2]; CHECK THIS AGAIN (*)
    P = [3e-6; rad^2];
    P = diag(P);
    
    %% Measurement matrix
    H = [-1 0];    
    
    %% measurment noise covariace matrix
    R = mag_noise_SD^2;
    
    for i = 1:size(time,1)
        %% Propogate state estimates
        x = Psi * x;
        %% Propogate error covariance matrix
        %P = P * [(scale_factor_error*DR_data(i))^2 0; 0 1]; CHECK THIS AGAIN (*)
        P = Psi * P * Psi' + Q;
        %% Compute Kalman Gain matrix
        K = P * H' * inv(H*P*H' + R);
        
        %% Measurment innovation vector
        z = compass(i) - gyro(i);
        dz = z - H*x;
        
        %% update state estimate
        x = x + K * dz;
        
        %% update error covariance matrix
        P = (eye(2) - K*H) * P;
        
        %% use state error estimate to correct heading solution
        heading(i) = gyro(i) - x(1);
    end
end