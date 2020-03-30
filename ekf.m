close all;
%clear;

% Model Parameters
T = 0.05;
C_dp = 1; % Drag coeff of parachute
A_p = 0; % Cross-sectional area of parachute (m^2)
C_dm = 0.25; % Drag coeff of package
A_m = 0.1016*0.0762; % Cross-sectional area of package (m^2)
m = 0.265; % Mass of package (assumes parachute is negligible)

% Read measurement file
data = load('data/drop2_long_long.mat').data;

% Kalman Parameters
Q = diag([0.01,0.01,0.1]);
R = diag([0.3814, 8.6526]);
p_kk = diag([R(1,1),0.01,0.01]);

% Logging variables
counter = 2;%12105;%round(356.4/0.05+1.5/0.05);
t = [0];
meas = [data(1, 2); data(1,5)-9.81];
acc_mag_log = [sqrt(data(1,3)^2+data(1,4)^2+data(1,5)^2)];

% State of system (which model to use)
parachute_time = 0; % Time in second at which the parachute opens
impact_time = 0;
impact_threshold = 50;
impact = false;
parachute_deployed = true;
free_fall_time = 0;
free_fall = false;

% Wait until free fall starts
while ~free_fall && counter <= length(data(:,1))
    % Check for freefall
    [free_fall, acc_mag] = is_free_fall(data(counter, :), data(counter-1,:), free_fall);
    
    % Log data
    z = [data(counter, 2); data(counter,5)-9.81];
    counter = counter +1;
    t = [t t(end)+T];
    meas = [meas z];
    acc_mag_log = [acc_mag_log acc_mag];
end

% Set free_fall_time, and initial ekf initial state
free_fall_time = t(end);
t_est = [t(end)];
x_kk = [z(1); 0; -9.81]; % m m/s m/s^2
x_est = [x_kk];

% Start estimator now that free fall has begun
while counter <= length(data(:,1))%12105+31%round(356.4/0.05 + 4/0.05)-35
    % Read measurements
    z = [data(counter, 2); data(counter,5)-9.81];
    
    % Check for parachute deployment
    if A_p == 0 && parachute_deployed &&  z(2)>=0
        parachute_time = t(end)+T;
        A_p = pi*(0.381/2)^2; % Open parachute
    end
   
    % Check for impact
    if ~impact
        [impact,acc_mag] = is_impact(data(counter,:), impact_threshold);
        if impact
            impact_time = t(end)+T;
            R = diag([0.3814, 1000]);
        end
    else
        % Log acc_mag data
        [~,acc_mag] = is_impact(data(counter,:), impact_threshold);
    end
    
    % Kalman Filter
    x_kkm1 = nonLinearModel(x_kk, T, m, C_dp, A_p, C_dm, A_m, impact);
    F = linearizedModel(x_kkm1, T, m, C_dp, A_p, C_dm, A_m, impact);
    p_kkm1 = F*p_kk*F' + Q;
    H = [1 0 0; 0 0 1];
    K = p_kkm1*H'*(H*p_kkm1*H'+R)^-1;
    x_kk = x_kkm1 + K*(z - H*x_kkm1);
    p_kk = (eye(3) - K*H)*p_kkm1;

    % Log data
    x_est = [x_est x_kk];
    t = [t t(end)+T];
    t_est = [t_est t_est(end)+T];
    meas = [meas z];
    acc_mag_log = [acc_mag_log acc_mag];
    counter = counter+1;
end

figure;
p1 = subplot(3,1,1);
plot(t_est, x_est(1,:), t, meas(1,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("Estimated height", "Height measurements");
ylabel("Height (m)");
grid on;

p2 = subplot(3,1,2);
plot(t_est, x_est(2,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("Estimated velocity");
ylabel("Velocity (m/s)");
grid on;

p3 = subplot(3,1,3);
plot(t_est, x_est(3,:), t, meas(2,:), t, acc_mag_log);
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
xline(impact_time,':', "Impact");
legend("Estimated acc", "Z acc measurements", "Acc magnitude measurements");
ylabel("Acceleration (m/s^2)");
grid on;
ylim([-50 150]);
linkaxes([p1,p2,p3],'x');
sgtitle("Extended Kalman Filter Parachute Drop Test");
xlabel('Time (s)');
export_fig -r500 'drop2_ekf.png'