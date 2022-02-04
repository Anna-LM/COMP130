function [range2sat, C_Ie_set] = compute_user2sat_ranges(cartECEF, cart_pos)
    % r_aj ranges from the approximate user position to each satellite
    % r_ej is cartesian ECEF position of the satellite j
    % r_ea is predicted Cartesian ECEF user position
    % C_Ie is the Sagnac effect compensation matrix

    r_ea = cartECEF;
    satNum = size(cart_pos, 2);
    
    C_Ie_set = zeros(3,3,satNum);
    range2sat = zeros([1,satNum]);
    
    C_Ie = eye(3);
    for j = 1:satNum
        r_ej = cart_pos(:,j);
        r_aj = sqrt(transpose(C_Ie*r_ej - r_ea) * (C_Ie*r_ej - r_ea));
        range2sat(1,j) = r_aj;
    end
    
    for j = 1:satNum
        C_Ie = sagnac_mat(range2sat(1,j));
        C_Ie_set(:,:,j) = C_Ie;
        r_ej = cart_pos(:,j);
        r_aj = sqrt(transpose(C_Ie * r_ej - r_ea) * (C_Ie * r_ej - r_ea));
        range2sat(1,j) = r_aj;
    end
end