function R = compute_measurement_noise_covmat(sat, rangeSD, rangerateSD)
    satNum = size(sat,2)-1;
    SDpr = rangeSD^2 * eye(satNum);
    SDprr = rangerateSD^2 * eye(satNum);
    Z = zeros(satNum);
    
    R = [SDpr Z;
         Z SDprr];
end