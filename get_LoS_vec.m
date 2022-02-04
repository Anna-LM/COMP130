function u_aj = get_LoS_vec(r_ea, r_ej, r_aj, C_Ie)
    u_aj = (C_Ie*r_ej - r_ea) / r_aj;
end