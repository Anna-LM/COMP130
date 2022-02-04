function u_aj = line_of_sight_vec(r_ea, r_aj, r_ej, C_Ie)
    
    num = C_Ie*r_ej - r_ea;

    u_aj = (num) / r_aj;


end