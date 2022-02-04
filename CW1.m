clear all;
clc;


% GNSS SOLUTION Only - data to use
pseudo_ranges = csvread('Pseudo_ranges.csv');
pseudo_rr = csvread('Pseudo_range_rates.csv');

[cartesian_position, innovation_vectors, measurement_matrices, range_from_p_to_sat, Line_of_sight_init, init_sagnac, clock_error] = GET_CART_ECEF_POS(pseudo_ranges);


Define_Constants;
outlier_present = 10;
while outlier_present > 0
    %Outlier detection for sateillites
    [pseudo_ranges, pseudo_rr, outlier_present] = Outlier_detection_SAT(pseudo_ranges, pseudo_rr, innovation_vectors, measurement_matrices);
    %Find position using pseudo ranges for lawnmower
    [cartesian_position, innovation_vectors, measurement_matrices, range_from_p_to_sat, Line_of_sight_init, init_sagnac, clock_error] = GET_CART_ECEF_POS(pseudo_ranges);
    
end
%Once Outliers are excluded pseudo range rates are calculated and in turn
%the 
[velocity_over_time, time_stamps, sattouser_range_r, clock_drifting] = GET_CART_ECEF_VEL(pseudo_ranges, pseudo_rr);


[GNSS_results, GNSS_display_table] = GET_NED_POSITION_VELOCITY(...
    cartesian_position, velocity_over_time, time_stamps);

%Altered csv code for GNSS_Only solution
dlmwrite('Initial_GNSS_solution.csv', GNSS_results, 'delimiter', ',', 'precision', 10.5);



% DEAD RECKONING SOLUTION
DR_data = csvread('Dead_reckoning.csv');

L0 = GNSS_results(1,2)*rad;
lamda0 = GNSS_results(1,3)*rad;
h0 = GNSS_results(1,4);

[DR_results, DR_display_table, avg_wheel_speed] = GET_DEAD_RECKONING_DATA(DR_data, L0, lamda0, h0);
dlmwrite('DR_only_solution.csv', DR_results, 'delimiter', ',', 'precision', 15);


% GNSS SOLUTION (with KALMAN FILTERING)
init_r_ea = transpose(cartesian_position(1,:));
init_v_ea = transpose(velocity_over_time(1,:));
init_clk_ofs = clock_error(1);
init_clk_dft = clock_drifting(1);
SD_ofs = 100000;
SD_dft = 200;
[x,P] = Init_GNSS_KF(init_r_ea, init_v_ea, init_clk_ofs, init_clk_dft, SD_ofs, SD_dft); 

tau = 0.5;
Phi = compute_GNSS_transition_matrix(tau);

accPSD = 0.01;
clkphsPSD = 0.01;
clkfreqPSD = 0.04;
Q = compute_GNSS_noise_covmat(tau, accPSD, clkphsPSD, clkfreqPSD);

rangeSD = 10;
rangerateSD = 0.05;

[GNSS_KF_results, GNSS_KF_display_table] = ...
    GET_KF_GNSS(pseudo_ranges, pseudo_rr, x, P, Phi, Q, rangeSD, rangerateSD);
dlmwrite('GNSSKF_solution.csv', GNSS_KF_results, 'delimiter', ',', 'precision', 15);



% INTEGRATED KF SOLUTION
mag_noise = 4*rad;          % magnetic compass noise SD
gyro_meas_noise = 1e-4;     % gyro noise SD
gyro_noise = 3e-6;          % gyro random noise SD
LC_KF_config.gyro_bias_PSD = 3e-6; % is gyro random walk noise
sf_err_SD = 0.01;           % scale factor error SD
wheel_speed_PSD = 0.01;


gyro = DR_data(:,6);
compass = DR_data(:,7)*deg_to_rad;

gyro_heading = compass(1);
result = zeros(size(gyro));

for i = 1:size(gyro,1)
    gyro_heading = gyro_heading + gyro(i) * tau;
    result(i) = gyro_heading;
end

gyro_heading_result = result;


final_results = GET_INTEGRATED_KF_SOLUTION(DR_results, gyro_heading_result, DR_data, GNSS_KF_results, tau, ...
    wheel_speed_PSD, gyro_meas_noise, gyro_noise, LC_KF_config.gyro_bias_PSD, mag_noise);
final_header = {'TIME(s)', 'LATITUDE(°)', 'LONGITUDE(°)', 'NORTH(m/s)', 'EAST(m/s)', 'HEADING(°)'};
final_display = [final_header; num2cell(final_results)];
dlmwrite('integrated_KF_solution.csv', final_results, 'delimiter', ',', 'precision', 15);

%% PLOTS
x0 = 50;
y0 = 50;
width = 500;
height = 750;

%% GNSS only plots
figure;
subplot(2,1,1);
plot(GNSS_results(:,3),GNSS_results(:,2));
title('NED Position of Lawnmower - Unfiltered GNSS');
xlabel('Longitude(°)');
ylabel('Latitude(°)');

subplot(2,1,2);
plot(GNSS_results(:,1),GNSS_results(:,5),GNSS_results(:,1),GNSS_results(:,6));
title('NED Velocity of Lawnmower - Unfiltered GNSS');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend({'North velocity','East velocity'}, 'Location','southeast');

set(gcf,'position',[x0,y0,width,height]);

%% KF GNSS plots
figure;
subplot(2,1,1);
plot(GNSS_KF_results(:,3),GNSS_KF_results(:,2));
title('NED Position of Lawnmower - Kalman Filtered GNSS');
xlabel('Longitude(°)');
ylabel('Latitude(°)');

subplot(2,1,2);
plot(GNSS_KF_results(:,1),GNSS_KF_results(:,5),GNSS_KF_results(:,1),GNSS_KF_results(:,6));
title('NED Velocity of Lawnmower - Kalman Filtered GNSS');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend({'North velocity','East velocity'}, 'Location','southeast');

set(gcf,'position',[x0,y0,width,height]);

%% DR plots
figure;
subplot(2,1,1);
plot(DR_results(:,3),DR_results(:,2));
title('NED Position of Lawnmower - Dead Reckoned');
xlabel('Longitude(°)');
ylabel('Latitude(°)');

subplot(2,1,2);
plot(DR_results(:,1),DR_results(:,4),DR_results(:,1),DR_results(:,5));
title('NED Velocity of Lawnmower - Dead Reckoned');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend({'North velocity','East velocity'}, 'Location','southeast');

set(gcf,'position',[x0,y0,width,height]);

%% Integrated Kalman Filter plots
figure;
subplot(3,1,1);
plot(final_results(:,3),final_results(:,2));
title('NED Position of Lawnmower - Integrated Kalman Filter');
xlabel('Longitude(°)');
ylabel('Latitude(°)');

subplot(3,1,2);
plot(final_results(:,1),final_results(:,4),final_results(:,1),final_results(:,5));
title('NED Velocity of Lawnmower - Integrated Kalman Filter');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend({'North velocity','East velocity'}, 'Location','southeast');

subplot(3,1,3);
plot(final_results(:,1),final_results(:,6));
title('Heading of Lawnmower - Integrated Kalman Filter');
xlabel('Time(s)');
ylabel('Heading(°)');

%% Comparison
figure;
plot(DR_results(:,3),DR_results(:,2),GNSS_KF_results(:,3),GNSS_KF_results(:,2),final_results(:,3),final_results(:,2));
title('NED Position of Lawnmower');
xlabel('Longitude(°)');
ylabel('Latitude(°)');
legend('Dead Reckoned solution','Kalman Filter GNSS Solution','Integrated Kalman Filter Solution');