function R = integrated_measurement_noise_covmatrix(SDGr, SDGv, L, h, mag_noiseSD)
    rad = pi/180;

    [RN,RE] = Radii_of_curvature(L);

    R = [SDGr^2/(RN+h)^2;
         SDGr^2/((RE+h)^2*cos(L)^2);
         SDGv^2;
         SDGv^2;
         mag_noiseSD^2;
         rad^2];
    
    R = diag(R);
end