function [in_ff, a_mag_c] = is_free_fall(current, prev, prev_ff)

    %% find magnitude of acceleration
    ff_mag =@(data) sqrt(data(3)^2+data(4)^2+data(5)^2);
    a_mag_c = ff_mag(current);
    a_mag_p = ff_mag(prev);
    
    %% logic
    low_thresh = 2;
    high_thresh = 4;
    in_ff = prev_ff;
    if  a_mag_c < low_thresh && a_mag_p < low_thresh
       in_ff = true;
    end
    if a_mag_c > high_thresh
        in_ff = false;
    end
end
