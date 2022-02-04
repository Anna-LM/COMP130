function r = get_user2sat_rangerate(r_ea, r_ej, u_aj, C_Ie, v_ea, v_ej)
    Define_Constants;
    r = u_aj' * (C_Ie *(v_ej + Omega_ie*r_ej) - (v_ea + Omega_ie*r_ea));
end