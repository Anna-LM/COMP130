function [cartesian_position, innovation_vectors, measurement_matrices, range_from_p_to_sat, Line_of_sight_init, init_sagnac, clock_error] = GET_CART_ECEF_POS(range_data, tolerance, position_offset, r_ea)
    %ECEF POSITIONS

    %Range data initialisation

    time_stamps = range_data(2:size(range_data,1),1);
    satellite = transpose(range_data(1,2:size(range_data,2)));
    cartesian_position = zeros(size(time_stamps,1), 3);
    clock_error = zeros(size(time_stamps,1), 1);
    range_from_p_to_sat = zeros(size(time_stamps,1), size(satellite,1));
    init_sagnac = zeros(3,3,size(satellite,1),size(time_stamps,1));
    Line_of_sight_init = zeros(3,size(satellite,1),size(time_stamps,1));
    data = range_data(2:size(range_data),2:size(range_data,2));
    innovation_vectors = zeros(size(satellite,1),size(time_stamps,1));
    measurement_matrices = zeros(size(satellite,1),4,size(time_stamps,1));
    origin_plane = [0;0;0];
    no_offset = 0;
    
    % Workshop 1 q1
    %determine the number of input arguments to position function in order
    %to determine set offset tolerance or range in cartesian plane
    
    if nargin <= 2
        
        r_ea = origin_plane;
        position_offset = no_offset;
        tolerance = 1e-6;
    elseif nargin <= 3
        r_ea = origin_plane;
        position_offset = no_offset;
    elseif nargin <= 4
        r_ea = origin_plane;
    elseif nargin > 5
        
        disp('error')
    end
        disp('')

    % Part b
    for i = 1:size(time_stamps,1)
        % arbitrary tolerance and difference
        difference = 100;
        while difference >= tolerance
            
            initialize_plane = zeros(3,size(satellite,1));
            sat_cartesian_plane = initialize_plane;

            for j = 1:size(satellite,1)
                %find the time and position of sateillite using sateillite
                %positiona and velocity function.

                r_ej = transpose(Satellite_position_and_velocity(time_stamps(i), satellite(j)));
                sat_cartesian_plane(:,j) = r_ej;

            % Part c
                [r_aj, C_Ie] = get_range_eq(r_ea, r_ej);
                range_from_p_to_sat(i,j) = r_aj;
                init_sagnac(:,:,j,i) = C_Ie;
                %initilaise sagnac equation and iterate for times

            % Part d
            %Compute the line of sight vector from user position to each
            %sateillite
                Line_of_sight_init(:,j,i) = line_of_sight_vec(r_ea, r_aj, r_ej, C_Ie);
            end

            % Part e
            x_hat = [r_ea; 
                position_offset];
            %normalize vector
            x_norm = norm(x_hat);
            dz = get_dz_vec(data(i,:), range_from_p_to_sat(i,:), position_offset);
%             H_g = [-1*Line_of_sight_init' ones(size(Line_of_sight_init),1)];
%             H_g = H_g(Line_of_sight_init(:,:,i));
            H_g = get_meas_matrix(Line_of_sight_init(:,:,i));

            %iterate for dz and H vectors
            innovation_vectors(:,i) = dz;
            measurement_matrices(:,:,i) = H_g;

            % Part f
            %update for time and position
            update = (inv(H_g'*H_g))*H_g'*dz;
            x_hat = x_hat + update;
            difference = abs(norm(x_hat) - x_norm);
            r_ea = x_hat(1:3);
            position_offset = x_hat(4);
            
        end
        %cartesian position and clock offset for times

        cartesian_position(i,:) = transpose(r_ea);
        clock_error(i) = position_offset;
    end
end