function [r, C_Ie] = get_range_eq(r_ea, r_ej)
    Define_Constants;
    
    r = sqrt(transpose(eye(3)*r_ej - r_ea) * (eye(3)*r_ej - r_ea));
    
    C_Ie = [1 omega_ie * r/c 0;
        -omega_ie * r/c 1 0;
        0 0 1];
    
    r = sqrt(transpose(C_Ie*r_ej - r_ea) * (C_Ie*r_ej - r_ea));
end