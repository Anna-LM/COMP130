function [results, display] = CLOSED_LOOP_INT_KF(start_lat, start_long, DR_data, GNSS_Kalman_f_data, tau_s, start_height, PSD)
    

    t = DR_data(:,1);
    comlumn_sum = sum(DR_data(:,4:5),2);

    GNSS_Lat = GNSS_Kalman_f_data(:,2)*rad;
    GNSS_lamda = GNSS_Kalman_f_data(:,3)*deg_to_rad;
    GNSS_h = GNSS_Kalman_f_data(:,4);
    GNSS_vN = GNSS_Kalman_f_data(:,5);
    GNSS_vE = GNSS_Kalman_f_data(:,6);

    mean_speed = comlumn_sum / 2;
    heading = DR_data(:,7)*deg_to_rad;
    
    lats = zeros(size(t));
    longitudes = lats;
    velocity_north = lats;
    velocity_east = lats;
    
      
    SDGv = 0.05; % wheel speed sensor noise
   
    % meas matrix
    

    H = [0 0 -1 0;
        0 0 0 -1;
        -1 0 0 0;
        0 -1 0 0];
    
    damped_constant = 1.7;
    damped_constant2 = 0.7;
    SDGr = 2;


    for i = 1:size(t,1)

        if i <= 1
            vN = (1/2) * (cos(heading(i)) + cos(heading(i)))*mean_speed(i);%avg velocity
            vE = (1/2) * (sin(heading(i)) + sin(heading(i)))*mean_speed(i);
       
        else
            vN = (1/2) * (cos(heading(i)) + cos(heading(i-1)))*mean_speed(i);
            vE = (1/2) * (sin(heading(i)) + sin(heading(i-1)))*mean_speed(i);


            disp('')

     
        end

        
        
        if i <= 1

            [RN, RE] = Radii_of_curvature(start_lat);

            L = start_lat + vN*(t(i)-t(i)) / (RN + start_height);
            lamda = start_long + vE*(t(i)-t(i)) / ((RE + start_height)*cos(L));
        else
            [RN, RE] = Radii_of_curvature(latitudes(i-1));
            L = latitudes(i-1) + (vN*(t(i)-t(i-1))) / (RN + start_height);
            lamda = longitudes(i-1) + vE*(t(i)-t(i-1)) / ((RE + start_height)*cos(L));
        end
        
        if i <= 1
            vN = avg_speed(1) * cos(heading(1));
            vE = avg_speed(1) * sin(heading(1));
        else
            vN = damped_constant*vN - damped_constant2*velocity_north(i-1);%damped instantaneous velocity - North and East
            vE = damped_constant*vE - damped_constant2*velocity_east(i-1);
        end
        
        x = [0;
            0;
            0;
            0];
        phi = compute_transition_mat(tau_s, L, GNSS_h(i));
        Q = compute_sys_noise_covmat(tau_s, PSD, L, GNSS_h(i));
        if i == 1
            % Initialise error covariance matrix
            [RN, RE] = Radii_of_curvature(L);
            P = [1; 1; 10^2/(RN+GNSS_h(1))^2; 10^2/((RE+GNSS_h(1))^2*cos(L))];
            P = diag(P);
        else
            P = phi * P * transpose(phi) + Q;
        end
        R = integrated_measurement_noise_covmat(SDGr, SDGv, L, GNSS_h(i));
        K = P * tranpose(H) * (H*P*transpose(H) + R);
        dz = [GNSS_Lat(i) - L;
             GNSS_lamda(i) - lamda;
             GNSS_vN(i) - vN;
             GNSS_vE(i) - vE];

        x = x + K*dz;
        
        L = L - x(3);
        lamda = lamda - x(4);
        vN = vN - x(1);
        vE = vE - x(2);
        
        latitudes(i) = L;
        longitudes(i) = lamda;
        velocity_north(i) = vN;
        velocity_east(i) = vE;
    end
    
    results = zeros(size(t,1),5);
    results(:,1) = t;
    results(:,2) = latitudes*deg;
    results(:,3) = longitudes*deg;
    results(:,4) = velocity_north;
    results(:,5) = velocity_east;
    
    header = {'time', 'KF latitiude', 'KF longitude', 'KF vel North', 'KF vel East'};
    display = [header; num2cell(results)];
end