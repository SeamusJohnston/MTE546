function F = linearizedModel(x_kkm1, T, m, C_dp, A_p, C_dm, A_m, impact)
%LINEARIZEDMODEL Linearized motion model matrix F
if ~impact
    F = [1 T 0.5*T^2; 0 1 T; 0 1.229*x_kkm1(2)/m*(C_dp*A_p + C_dm*A_m) 0];
else
    F = [1 T 0.5*T^2; 0 1 T; 0 -1/T 0];
end

