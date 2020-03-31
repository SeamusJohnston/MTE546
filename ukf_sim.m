close all;

% Model Parameters
T = 0.05;
C_dp = 1; % Drag coeff of parachute2
A_p = 0; % Cross-sectional area of parachute (m^2)
C_dm = 0.25; % Drag coeff of package
A_m = 0.1016*0.0762; % Cross-sectional area of package (m^2)
m = 0.265; % Mass of package (assumes parachute is negligible)

% Kalman Parameters
L = 3;
alpha=1e-3;
ki=0;
beta=2;
Q = diag([0.0001,0.0001,0.1]);
R = diag([0.3814, 0.1]);
p_kk = diag([0.3814,0.01,0.01]);

% Logging variables
counter = 1;
t = [0];
meas = [100; -9.81];

% State of system (which model to use)
parachute_time = 0; % Time in second at which the parachute opens
parachute_deployed = true;

% Set free_fall_time, and initial ekf initial state and truth state
free_fall_time = 0;
t_est = [t(end)];
x_kk = [100; 0; -9.81]; % m m/s m/s^2
x_est = [x_kk];
x = x_est;
x_truth = x_kk;

while counter <= 50
    if counter == 10
        parachute_time = t(end)+T;
        A_p = pi*(0.381/2)^2; % Open parachute
    end
    
    % Read measurements and compute truth data
    x_truth = nonLinearModel(x_truth, T, m, C_dp, A_p, C_dm, A_m, impact)+ [sqrt(Q(1,1)); sqrt(Q(2,2)); sqrt(Q(3,3))].*randn(3,1);
    z = [x_truth(1); x_truth(3)]+ [sqrt(R(1,1)); sqrt(R(2,2))].*randn(2,1);
    
    % Kalman Filter
    [X, Wm, Wc] = sigmaPoints(x_kk, p_kk, alpha, beta, ki);
    X_f = X;
    for i=1:numel(X(1,:))
        X_f(:,i) = nonLinearModel(X(:,i), T, m, C_dp, A_p, C_dm, A_m, impact);
    end
    x_kkm1 = sum(Wm.*X_f,2);
    p_kkm1 = Wc.*(X_f-x_kkm1)*(X_f-x_kkm1)' + Q;
    [X, Wm, Wc] = sigmaPoints(x_kk, p_kk, alpha, beta, ki);
    H = [1 0 0; 0 0 1];
    Z = zeros(2,numel(X(1,:)));
    for i=1:numel(X(1,:))
        Z(:,i) = H*X(:,i);
    end
    z_bar = sum(Wm.*Z,2);
    S = Wc.*(Z-z_bar)*(Z-z_bar)' + R;
    C_sz = Wc.*(X_f-x_kkm1)*(Z-z_bar)'; % Cross cov
    K = C_sz*(S)^-1;
    x_kk = x_kkm1 + K*(z - z_bar);
    p_kk = p_kkm1 - K*S*K';

    % Log data
    x_est = [x_est x_kk];
    x= [x x_truth];
    t = [t t(end)+T];
    t_est = [t_est t_est(end)+T];
    meas = [meas z];
    counter = counter+1;
end

% Plot estimate
figure('Renderer', 'painters', 'Position', [10 10 300 400])
p1 = subplot(3,1,1);
plot(t_est, x_est(1,:), t, x(1,:), t, meas(1,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
%xline(impact_time,':', "Impact");
legend("Estimated", "Truth", "Measured");
ylabel("Height (m)");
grid on;

p2 = subplot(3,1,2);
plot(t_est, x_est(2,:), t, x(2,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
%xline(impact_time,':', "Impact");
legend("Estimated", "Truth");
ylabel("Velocity (m/s)");
grid on;

p3 = subplot(3,1,3);
plot(t_est, x_est(3,:),t, x(3,:), t, meas(2,:));
xline(free_fall_time,':', "Free Fall");
xline(parachute_time,':', "Parachute open");
%xline(impact_time,':', "Impact");
legend("Estimated", "Truth", "Measured");
ylabel("Acceleration (m/s^2)");
grid on;
%ylim([-50 150]);
linkaxes([p1,p2,p3],'x');
sgtitle("Unscented Kalman Filter Simulation");
xlabel('Time (s)');
set(gcf,'Color',[1 1 1])
export_fig -r500 'sim_ukf.png'
