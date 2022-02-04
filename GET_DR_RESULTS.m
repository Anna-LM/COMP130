function [results, display_table, avg_speed] = GET_DR_RESULTS(data, init_latitude, init_longitude, init_height)
    
    Define_Constants;

    times = data(:,1);
    avg_speed = sum(data(:,4:5),2) / 2;
    gyro = data(:,6);
    heading = data(:,7)*rad;

    [vN, vE] = get_NED_speed(avg_speed, heading);
    [DR_lat, DR_long] = get_DR_latitudes_longitudes...
        (init_latitude, init_longitude, init_height, vN, vE, times);

    v0 = avg_speed(1);
    Psi0 = heading(1);
    [vN, vE] = get_damped_instant_velocity(vN, vE, v0, Psi0);
    
    results = zeros(size(times,1),5);
    results(:,1) = times;
    results(:,2) = DR_lat*deg;
    results(:,3) = DR_long*deg;
    results(:,4) = vN;
    results(:,5) = vE;
    
    header = {'TIME(s)', 'DR LATITUDE(°)', 'DR LONGITUDE(°)', 'DR NORTH(m/s)', 'DR EAST(m/s)'};
    display_table = [header; num2cell(results)];
end