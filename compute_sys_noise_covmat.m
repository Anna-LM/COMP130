function Q = compute_sys_noise_covmat(tau, S, L, h)
    [RN,RE] = Radii_of_curvature(L);

    Q = [S*tau;
         S*tau; 
         ((1/3)*S*tau^3)/(RN+h)^2;
         ((1/3)*S*tau^3)/((RE+h)^2*cos(L)^2)];
     Q = diag(Q);
     Q(1,3) = (0.5*S*tau^2)/(RN+h);
     Q(2,4) = (0.5*S*tau^2)/((RE+h)*cos(L));
     Q(3,1) = Q(1,3);
     Q(4,2) = Q(2,4);
end