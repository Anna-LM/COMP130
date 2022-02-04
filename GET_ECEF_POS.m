function [satuser_range, cartesian_results, innovation_vectors, measurement_matrices, LoS_set, sagnac_set, clk_offsets] = ...
    GET_ECEF_POS(range_data, tolerance, offset, r_ea)

    t = range_data(2:size(range_data,1),1);
    sat = transpose(range_data(1,2:size(range_data,2)));
    
    
    satuser_range = zeros(size(t,1), size(sat,1));
    sagnac_set = zeros(3,3,size(sat,1),size(t,1));
    LoS_set = zeros(3,size(sat,1),size(t,1));
    cartesian_results = zeros(size(t,1), 3);
    clk_offsets = zeros(size(t,1), 1);
    
    data = range_data(2:size(range_data),2:size(range_data,2));
    innovation_vectors = zeros(size(sat,1),size(t,1));
    measurement_matrices = zeros(size(sat,1),4,size(t,1));

    % 1a.a)
    if nargin < 2
        r_ea = [0;0;0];
        offset = 0;
        tolerance = 1e-8;
    elseif nargin < 3
        r_ea = [0;0;0];
        offset = 0;
    elseif nargin < 4
        r_ea = [0;0;0];
    end

    %% 1a.b)
    for i = 1:size(t,1)
        %% reset difference each time
        % (1000 is a arbitrary number just to ensure difference is higher in
        % the beginning)
        diff = 1000;
        while diff > tolerance
            sat_cart = zeros(3,size(sat,1));
            for j = 1:size(sat,1)
                r_ej = transpose(Satellite_position_and_velocity(t(i), sat(j)));
                sat_cart(:,j) = r_ej;

            %% 1a.c)
                [r_aj, C_Ie] = get_user2sat_range(r_ea, r_ej);
                satuser_range(i,j) = r_aj;
                sagnac_set(:,:,j,i) = C_Ie;
            %% 1a.d)
                LoS_set(:,j,i) = get_LoS_vec(r_ea, r_ej, r_aj, C_Ie);
            end

            %% 1a.e)
            x = [r_ea; offset];
            x_norm = norm(x);
            dz = get_meas_innovation_vec(data(i,:), satuser_range(i,:), offset);
            H = get_meas_mat(LoS_set(:,:,i));

            %% store innovation vector and measurement matrix
            innovation_vectors(:,i) = dz;
            measurement_matrices(:,:,i) = H;

            %% 1a.f)
            x = calculate_new_state(x, H, dz);
            r_ea = x(1:3);
            offset = x(4);
            diff = abs(norm(x) - x_norm);
        end

        cartesian_results(i,:) = transpose(r_ea);
        clk_offsets(i) = offset;
    end
end