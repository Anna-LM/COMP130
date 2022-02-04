function [results, display_table] = GET_NED_POSITION_VELOCITY(cartesian_positions, cartesian_velocities, times)
    Define_Constants;
    results = zeros(size(times,1), 7);
    for i = 1:size(times, 1)
        [L, lamda, h, vel] = pv_ECEF_to_NED(transpose(cartesian_positions(i,:)), transpose(cartesian_velocities(i,:)));
        results(i,1) = times(i);
        results(i,2) = L*deg;
        results(i,3) = lamda*deg;
        results(i,4) = h;
        results(i,5:7) = vel';
    end
    
    header = {'TIME(s)', 'LATITUDE(°)', 'LONGITUDE(°)', 'HEIGHT(m)', 'NORTH(m/s)', 'EAST(m/s)', 'DOWN(m/s)'};
    display_table = [header; num2cell(results)];
end