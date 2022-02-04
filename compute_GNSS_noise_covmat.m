function Q = compute_GNSS_noise_covmat(tau, accPSD, clkphsPSD, clkfreqPSD)
    Q = zeros(8);
    Q(1:3,1:3) = (1/3) * accPSD * tau^3 * eye(3);
    Q(1:3,4:6) = 0.5 * accPSD * tau^2 * eye(3);
    Q(4:6,1:3) = Q(1:3,4:6);
    Q(4:6,4:6) = accPSD * tau * eye(3);
    Q(7,7) = clkphsPSD * tau + (1/3) * clkfreqPSD * tau^3;
    Q(7,8) = 0.5 * clkfreqPSD * tau^2;
    Q(8,7) = Q(7,8);
    Q(8,8)= clkfreqPSD * tau;
end