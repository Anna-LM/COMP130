function Phi = compute_GNSS_transition_matrix(propogation_interval)
    tau = propogation_interval;
    
    Phi = eye(8);
    Phi(1:3,4:6) = tau * eye(3);
    Phi(7,8) = tau;
end