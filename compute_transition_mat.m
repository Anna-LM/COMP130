function Psi = compute_transition_mat(tau, L, h)
    Psi = eye(4);
    [RN,RE] = Radii_of_curvature(L);
    
    Psi(3,1) = tau/(RN + h);
    Psi(4,2) = tau/((RE + h)*cos(L));
end