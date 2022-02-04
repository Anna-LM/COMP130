function [x, P] =  Init_GNSS_KF(init_r_ea, init_v_ea, init_clk_ofs, init_clk_dft, SD_ofs, SD_dft)
    %% Initialise Kalman Filter state vector estimate
    x = [init_r_ea;
         init_v_ea;
         init_clk_ofs;
         init_clk_dft];         

    %% Initialise error covariance matrix
    P = [SD_ofs^2*ones(3,1);
         SD_dft^2*ones(3,1);
         SD_ofs^2;
         SD_dft^2];
    P = diag(P);
end