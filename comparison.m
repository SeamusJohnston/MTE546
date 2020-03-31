% Superimpose plots
clear;
ekf;
x_est_ekf = x_est;
ukf;
x_est_ukf = x_est;

close all;

figure('Renderer', 'painters', 'Position', [10 10 300 400]);
p1 = subplot(3,1,1);
plot(t_est, x_est_ukf(1,:), t_est, x_est_ekf(1,:), height_gt(:,1)-0.3, height_gt(:,2), t, meas(1,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("UKF", "EKF", "Truth", "Measured");
ylabel("Height (m)");
grid on;

p2 = subplot(3,1,2);
plot(t_est, x_est_ukf(2,:), t_est, x_est_ekf(2,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("UKF","EKF");
ylabel("Velocity (m/s)");
grid on;

p3 = subplot(3,1,3);
plot(t_est, x_est_ukf(3,:), t_est, x_est_ekf(3,:), t, meas(2,:), t, acc_mag_log);
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("UKF", "EKF", "Measured", "Magnitude");
ylabel("Acceleration (m/s^2)");
grid on;
ylim([-50 150]);
linkaxes([p1,p2,p3],'x');
sgtitle("Summary Kalman Filter Drop Test");
xlabel('Time (s)');
set(gcf,'Color',[1 1 1])
export_fig -r500 'drop2_comp.png'

% Plot the error between ground truth and the estimates
ukf_height = x_est_ukf(1,1:27);
ekf_height = x_est_ekf(1,1:27);
t_est = t_est(1:27);
height = [height_gt(10:end-2,1)-0.3 height_gt(10:end-2,2)];
t_i = 0.45:0.1:1.75;

true_height_i = interp1(height(:,1), height(:,2), t_i);
ukf_height_i = interp1(t_est, ukf_height, t_i);
ekf_height_i = interp1(t_est, ekf_height, t_i);

figure;
plot(t_i, abs(ukf_height_i-true_height_i))
hold on;
plot(t_i, abs(ekf_height_i-true_height_i))
grid on;
ylabel("Absolute Error (m)")
xlabel("Time(s)")
title("Estimation Error in Altitude")
ylim([0 0.45]);
legend("UKF", "EKF");