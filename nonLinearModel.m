function x_kkm1 = nonLinearModel(x_km1km1, T, m, C_dp, A_p, C_dm, A_m)
%nonLinearModel Non linear motion model
x_kkm1 = zeros(3,1);
x_kkm1(1) = x_km1km1(1) + x_km1km1(2)*T + 0.5*x_km1km1(3)*T^2;
x_kkm1(2) = x_km1km1(2) + x_km1km1(3)*T;
x_kkm1(3) = 1.229*x_km1km1(2)^2/(2*m)*(C_dp*A_p + C_dm*A_m) - 9.81;
end

