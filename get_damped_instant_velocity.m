function [vN, vE] = get_damped_instant_velocity(vN_, vE_, v0, heading0)
    vN = zeros(size(vN_));
    vE = vN;

    for i = 1:size(vN_,1)
        if i == 1
            vN(i) = v0 * cos(heading0);
            vE(i) = v0 * sin(heading0);
        else
            vN(i) = 1.7*vN_(i) - 0.7*vN(i-1);
            vE(i) = 1.7*vE_(i) - 0.7*vE(i-1);
        end
    end
end