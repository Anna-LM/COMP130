function [results, display, mean_speed] = GET_DEAD_RECKONING_DATA...
    (Get_dr_data, init_latitude, init_longitude, init_height)
    
    Define_Constants;

    times = Get_dr_data(:,1);
    mean_speed = sum(Get_dr_data(:,4:5),2) / 2;
    heading = Get_dr_data(:,7)*deg_to_rad;
    
    vN = zeros(size(mean_speed));
    vE = vN;
    %average velocity task 1 - Workshop 3
    %using equations exactly for mean speed
    for i = 1:size(mean_speed,1)
        if i == 1
            mean = 1/2;
            %initial speed
            vN(i) = mean * (cos(heading(i)) + cos(heading(i))) * mean_speed(i);
            vE(i) = mean * (sin(heading(i)) + sin(heading(i))) * mean_speed(i);

        elseif i> 1
            %velocity for front and back wheels
            vN(i) = mean * (cos(heading(i)) + cos(heading(i-1))) * mean_speed(i);
            vE(i) = mean * (sin(heading(i)) + sin(heading(i-1))) * mean_speed(i);
            
        elseif i < 1
            %error catcher most likely unnesscessary
            disp('error')

        end
    end
    
    %must find latitude and longitude on plane
    latitudes = zeros(size(times));
    longitudes = latitudes;
    %latititude and longitude equations from task 1 - Workshop 3
    for i = 1:size(times,1)
        if i == 1
            [RN, RE] = Radii_of_curvature(init_latitude);
            latitudes(i) = init_latitude + vN(i)*(times(i)-times(i)) / (RN + init_height);
            longitudes(i) = init_longitude + vE(i)*(times(i)-times(i)) / ((RE + init_height)*cos(latitudes(i)));
        elseif i>1
            [RN, RE] = Radii_of_curvature(latitudes(i-1));
            latitudes(i) = latitudes(i-1) + vN(i)*(times(i)-times(i-1)) / (RN + init_height);
            longitudes(i) = longitudes(i-1) + vE(i)*(times(i)-times(i-1)) / ((RE + init_height)*cos(latitudes(i)));
        elseif i < 1
            disp('error')
        end
    end

    %task 1 - final part Workshop 3
    %Get damped instantaneous DR velocity
    
    v0 = mean_speed(1);
    Psi0 = heading(1);
    [vN, vE] = get_damped_instant_velocity(vN, vE, v0, Psi0);
    
    results = zeros(size(times,1),5);
    results(:,1) = times;
    results(:,2) = latitudes*rad_to_deg;
    results(:,3) = longitudes*rad_to_deg;
    results(:,4) = vN;
    results(:,5) = vE;
    
    header = {'TIME(s)', 'DR LATITUDE(°)', 'DR LONGITUDE(°)', 'DR NORTH(m/s)', 'DR EAST(m/s)'};
    display = [header; num2cell(results)];
end