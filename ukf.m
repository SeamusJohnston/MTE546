close all;
clear;

% Model Parameters
T = 0.05;
C_dp = 1.75; % Drag coeff of parachute
A_p = 0; % Cross-sectional area of parachute
C_dm = 1; % Drag coeff of package
A_m = 0.1; % Cross-sectional area of package
m = 1; % Mass of package (assumes parachute is negligible)

% Kalman Parameters
L = 3;
alpha=1e-3;
ki=0;
beta=2;
Q = diag([0.01,0.01,0.01]);
R = diag([0.1, 0.5]);
x_kk = [500; 0; -9.81]; % m m/s m/s^2
p_kk = diag([0.01,0.01,0.01]);
x_real = x_kk;

% Logging variables
counter = 1;
t = [0];
x_est = [x_kk];
x_truth = [x_real];
meas = [x_kk(1); x_kk(3)];
parachute_time = 0; % Time in second at which the parachute opens

while counter <= 200
    if counter == 50
        parachute_time = t(end)+T;
        A_p = 0.1; % Open parachute
    end

    % Run Sim
    x_real = nonLinearModel(x_real, T, m, C_dp, A_p, C_dm, A_m);
    z = [x_real(1); x_real(3)] + [sqrt(R(1,1))*randn(1); sqrt(R(2,2))*randn(1)]; % Add noise
        
    % Kalman Filter
    [X, Wm, Wc] = sigmaPoints(x_kk, p_kk, alpha, beta, ki);
    X_f = X;
    for i=1:numel(X(1,:))
        X_f(:,i) = nonLinearModel(X(:,i), T, m, C_dp, A_p, C_dm, A_m);
    end
    x_kkm1 = sum(Wm.*X_f,2);
    p_kkm1 = Wc.*(X_f-x_kkm1)*(X_f-x_kkm1)' + Q;
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
    t = [t t(end)+T];
    x_truth = [x_truth x_real];
    meas = [meas z];
    counter = counter+1;
end

figure;
p1 = subplot(3,1,1);
plot(t, x_est(1,:), t, x_truth(1,:), t, meas(1,:));
xline(parachute_time,':', "Parachute open");
legend("Estimated height", "True height", "Height measurements");
ylabel("Height (m)");

p2 = subplot(3,1,2);
plot(t, x_est(2,:), t, x_truth(2,:));
xline(parachute_time,':', "Parachute open");
legend("Estimated velocity", "True velocity");
ylabel("Velocity (m/s)");

p3 = subplot(3,1,3);
plot(t, x_est(3,:), t, x_truth(3,:), t, meas(2,:));
xline(parachute_time,':', "Parachute open");
legend("Estimated acceleration", "True acceleration", "Acceleration measurements");
ylabel("Acceleration (m/s^2)");
linkaxes([p1,p2,p3],'x');
sgtitle("Unscented Kalman Filter Parachute Simulation");
xlabel('Time (s)');