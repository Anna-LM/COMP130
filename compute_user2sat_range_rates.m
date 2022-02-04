function rdot = compute_user2sat_range_rates(cartECEF, cartECEFvel, cart_pos, cart_vel, C_Ie_set, LoS_vec)
    Define_Constants;
    r_ea = cartECEF;
    v_ea = cartECEFvel;
    
    satNum = size(cart_pos, 2);
    rdot = zeros([1, satNum]);
    
    for i = 1:satNum
        r_ej = cart_pos(:,i);        % relative position of satellite j
        v_ej = cart_vel(:,i);        % relative velocity of satellite j
        u_aj = LoS_vec(:,i);         % Line of sight vectors for time t
        C_Ie = C_Ie_set(:,:,i);      % Sagnac matrix set for each satellite at time t
        
        rdot_aj = u_aj' * (C_Ie * (v_ej + Omega_ie*r_ej) - (v_ea + Omega_ie*r_ea));
        rdot(:,i) = rdot_aj;        % satellite range rate
    end
end