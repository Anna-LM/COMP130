function [LoS, LoS_approx] = calc_line_of_sight_vec(cartECEF, cart_pos, r_aj_set, C_Ie_set)
    r_ea = cartECEF;
    satNum = size(r_aj_set, 2);
    
    LoS = zeros([3, satNum]); % line of sight set
    LoS_approx = zeros([3, satNum]);
    
    for i = 1:satNum
        r_aj = r_aj_set(:,i);
        r_ej = cart_pos(:,i);
        C_Ie = C_Ie_set(:,:,i);
        
        % actual
        u_aj = (C_Ie*r_ej - r_ea) / r_aj;
        LoS(:,i) = u_aj;
        
        % approximation
        u_aj = (r_ej - r_ea) / r_aj;
        LoS_approx(:,i) = u_aj;
    end
end