% Superimpose plots
clear;
ekf;
x_est_ekf = x_est;
ukf;
x_est_ukf = x_est;

close all;

figure;
p1 = subplot(3,1,1);
plot(t_est, x_est_ukf(1,:), t_est, x_est_ekf(1,:), t, meas(1,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("Estimated height - UKF", "Estimated height - EKF", "Height measurements");
ylabel("Height (m)");
grid on;

p2 = subplot(3,1,2);
plot(t_est, x_est_ukf(2,:), t_est, x_est_ekf(2,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("Estimated velocity - UKF","Estimated velocity - EKF");
ylabel("Velocity (m/s)");
grid on;

p3 = subplot(3,1,3);
plot(t_est, x_est_ukf(3,:), t_est, x_est_ekf(3,:), t, meas(2,:), t, acc_mag_log);
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("Estimated acc - UKF", "Estimated acc - EKF", "Z acc measurements", "Acc magnitude measurements");
ylabel("Acceleration (m/s^2)");
grid on;
ylim([-50 150]);
linkaxes([p1,p2,p3],'x');
sgtitle("Summary Kalman Filter Parachute Test");
xlabel('Time (s)');
export_fig -r500 'drop2_comp.png'