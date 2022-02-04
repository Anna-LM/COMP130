function [vN, vE] = get_NED_speed(speed, heading)
    vN = zeros(size(speed));
    vE = vN;

    for i = 1:size(speed,1)
        if i == 1
            vN(i) = 0.5 * (cos(heading(i)) + cos(heading(i))) * speed(i);
            vE(i) = 0.5 * (sin(heading(i)) + sin(heading(i))) * speed(i);
        else
            vN(i) = 0.5 * (cos(heading(i)) + cos(heading(i-1))) * speed(i);
            vE(i) = 0.5 * (sin(heading(i)) + sin(heading(i-1))) * speed(i);
        end
    end
end