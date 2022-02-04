function Q = compute_sys_noise_covmatrix(tau, DR_PSD, L, h, gyro_noise, gyro_bias)
    [RN,RE] = Radii_of_curvature(L);

    Q = [DR_PSD*tau;
         DR_PSD*tau; 
         ((1/3)*DR_PSD*tau^3)/(RN+h)^2;
         ((1/3)*DR_PSD*tau^3)/((RE+h)^2*cos(L)^2);
         gyro_noise*tau+gyro_bias/3*tau^3;
         gyro_bias*tau];
     Q = diag(Q);
     Q(1,3) = (0.5*DR_PSD*tau^2)/(RN+h);
     Q(2,4) = (0.5*DR_PSD*tau^2)/((RE+h)*cos(L));
     Q(3,1) = Q(1,3);
     Q(4,2) = Q(2,4);
     Q(5,6) = 0.5*gyro_bias*tau^2;
     Q(6,5) = 0.5*gyro_bias*tau^2;
end