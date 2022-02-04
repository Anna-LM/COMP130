function C = sagnac_mat(r)
    omega = 7.292115e-5; %rad/s and is the Earth rotation rate
    c = 299792458; %m/s and is the speed of light
    
    %Sagnac effect compensation matrix
    C = [1 omega*r/c 0;
         -omega*r/c 1 0;
         0 0 1];
end