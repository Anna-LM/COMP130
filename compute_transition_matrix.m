function Psi = compute_transition_matrix(tau, L, h)
    Psi = eye(6);
    [RN,RE] = Radii_of_curvature(L);
    
    Psi(3,1) = tau/(RN + h);
    Psi(4,2) = tau/((RE + h)*cos(L));
    Psi(5,6) = tau;
end