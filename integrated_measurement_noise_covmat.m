function R = integrated_measurement_noise_covmat(SDGr, SDGv, L, h)
    [RN,RE] = Radii_of_curvature(L);

    R = [SDGr^2/(RN+h)^2;
         SDGr^2/((RE+h)^2*cos(L)^2);
         SDGv^2;
         SDGv^2];
    
    R = diag(R);
end