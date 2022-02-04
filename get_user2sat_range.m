function [r, C_Ie] = get_user2sat_range(user_cart, sat_cart)
    r_ea = user_cart;
    r_ej = sat_cart;

    r = sqrt(transpose(eye(3)*r_ej - r_ea) * (eye(3)*r_ej - r_ea));
    C_Ie = get_sagnac_mat(r);
    r = sqrt(transpose(C_Ie*r_ej - r_ea) * (C_Ie*r_ej - r_ea));
end