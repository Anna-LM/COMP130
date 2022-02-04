function [velocity_over_time, time_stamps, sattouser_range_r, clock_drifting] = GET_CART_ECEF_VEL(range_data, rangerate_data, v_ea, clock_drift)
    [cartesian_position, ~, measurement_matrices, ~, Line_of_sight_init, sagnac_set]= GET_CART_ECEF_POS(range_data);

    time_stamps = range_data(2:size(range_data,1),1);
    sateillite = transpose(range_data(1,2:size(range_data,2)));
    rangerate_data = rangerate_data(2:size(rangerate_data),2:size(rangerate_data,2));
    sattouser_range_r = zeros(size(time_stamps,1), size(sateillite,1));
    velocity_over_time = zeros(size(time_stamps,1), 3);
    clock_drifting = zeros(size(time_stamps,1), 1);
    set_clock_drift = 0;
    origin = [0;0;0];
    %Same approach as position
    %iterate through times and sateillites
    %do not need to include difference loop

    if nargin <= 3
        v_ea = origin;
        clock_drift = set_clock_drift;
    elseif nargin <= 4
        clock_drift = set_clock_drift;
    elseif nargin > 4
        disp('')


    end
    
    for i = 1:size(time_stamps,1)
        %initialise zero arrays for position and velocity 
        %in order to iterate through times
        pos_mat = zeros(3, size(sateillite,1));
        vel_mat = zeros(3, size(sateillite,1));

        %Use sateillite and velocity function 
        for j = 1:size(sateillite,1)
            [r_ej, v_ej] = Satellite_position_and_velocity(time_stamps(i), sateillite(j));
            r_ej = transpose(r_ej);
            v_ej = transpose(v_ej);
            %fill up zeros array

            r_ea = transpose(cartesian_position(i,:));

            pos_mat(:,j) = r_ej';
            vel_mat(:,j) = v_ej';
            
            %Line of sight and sagnac matrix
            LOS = Line_of_sight_init(:,j,i);
            SAG_MAT = sagnac_set(:,:,j,i);
            %predicted range rates equation

            Define_Constants;
            rdot_aj = transpose(LOS) * (SAG_MAT *(v_ej + Omega_ie*r_ej) - (v_ea + Omega_ie*r_ea));

            
            sattouser_range_r(i,j) = rdot_aj;        
        end
        %same approach as position - Workshop 1 
        x_hat = [v_ea; clock_drift];
        dz = get_dz_vec(rangerate_data(i,:), sattouser_range_r(i,:), clock_drift);
        H_g = measurement_matrices(:,:,i);
        
        %new state -update

        update = (inv(H_g'*H_g))*H_g'*dz;
        x_hat = x_hat + update;
        
        %velocity with clock drift
        v_ea = x_hat(1:3);
        clock_drift = x_hat(4);
    end

    %cartesian velocity and clock drift over time

    velocity_over_time(i,:) = transpose(v_ea);
    clock_drifting(i) = clock_drift;

end