function [results, display_table] = GET_KF_GNSS(range_data, rangerate_data, x, P, Phi, Q, rangeSD, rangerateSD)
    time = range_data(2:size(range_data,1), 1);
    results = zeros([size(time,1), 6]);

    for t = 1:size(time,1)
        [x, P, lat, long, height, vel] = ...
            Task2A_Function(range_data, rangerate_data, x, P, Phi, ...
                            Q, time(t), rangeSD, rangerateSD);
        results(t,1) = lat;
        results(t,2) = long;
        results(t,3) = height;
        results(t,4) = vel(1);
        results(t,5) = vel(2);
        results(t,6) = vel(3);
    end

    results = [time results];
    header = {'TIME(s)', 'LATITUDE(°)', 'LONGITUDE(°)', 'HEIGHT(m)', 'NORTH(m/s)', 'EAST(m/s)', 'DOWN(m/s)'};
    display_table = [header; num2cell(results)];
end