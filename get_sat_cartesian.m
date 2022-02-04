function [sat_cart_pos, sat_cart_vel] = get_sat_cartesian(data, t)
    satNums = data(1, 2:size(data,2));       %Satellite numbers
    sat_cart_pos = zeros([size(data,2)-1,3]);    %Cartesian postions vector
    sat_cart_vel = zeros([size(dataset,2)-1,3]); %Cartesian velocities vector
    
    for j = 1:size(satNums, 2)
        [sat_cart_pos(j,:), sat_cart_vel(j,:)] = Satellite_position_and_velocity(t, satNums(j));
    end
    
    sat_cart_pos = transpose(sat_cart_pos);
    sat_cart_vel = transpose(sat_cart_vel);
end