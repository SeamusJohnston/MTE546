function [impact, mag] = is_impact(accel_data,threshold)
%IS_IMPACT Determine if there has been an impact
impact = false;
mag = sqrt(accel_data(3)^2+accel_data(4)^2+accel_data(5)^2);
if mag >= threshold
    impact = true;
end

