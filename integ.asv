clc; clear; close all;

% Load UDDS data
load('C:\Users\USER\Desktop\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
udds_current = meas.Current; % UDDS 전류 데이터
udds_voltage = meas.Voltage; % UDDS 전압 데이터
udds_time = meas.Time; % UDDS 시간 데이터

% Load SOC-OCV structure
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);
ocv_values = soc_ocv(:, 2);

% Configuration parameters
Config.dt = mean(diff(udds_time)); % 평균 시간 간격
Config.ik = udds_current(1); % 초기 전류
Config.R0 = 0.001884314;
Config.R1 = 0.045801322;
Config.C1 = 4846.080679;
Config.cap = (2.99 * 3600) / 100; % nominal capacity [Ah]
Config.coulomb_efficient = 1;

% Initialize SOC estimation
SOC_est = zeros(length(udds_current), 1);
SOC_est(1) = 1; % 초기 SOC는 100%로 가정

% Initialize true SOC using coulomb counting
true_SOC = zeros(length(udds_current), 1);
true_SOC(1) = 1; % 초기 SOC는 100%로 가정

% EKF 초기화
P = eye(2); % Initial estimation error covariance

% Preallocate arrays for V1 and Vt
V1_est = zeros(length(udds_current), 1);
Vt_est = zeros(length(udds_current), 1);
V1_est(1) = udds_current(1) * Config.R1 * (1 - exp(-Config.dt / (Config.R1 * Config.C1))); % 초기 V1 값
Vt_est(1) = udds_voltage(1);

% Simulation loop
for k = 2:length(udds_current)
    Config.ik = udds_current(k);

    % True SOC calculation using coulomb counting
    delta_t = udds_time(k) - udds_time(k-1);
    true_SOC(k) = true_SOC(k-1) - (udds_current(k) * delta_t) / (Config.cap * 3600);

    % SOC estimation using EKF
    [SOC_est(k), V1_est(k), Vt_est(k), P] = soc_estimation(SOC_est(k-1), V1_est(k-1), udds_voltage(k), Config.ik, Config, P, soc_values, ocv_values);
end

% Plot SOC
figure;
plot(udds_time, true_SOC, 'b', 'LineWidth', 1.5); hold on;
plot(udds_time, SOC_est, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('SOC');
title('True SOC vs Estimated SOC during UDDS Cycle');
legend('True SOC', 'Estimated SOC');
grid on;

% Plot Vt
figure;
plot(udds_time, udds_voltage, 'b', 'LineWidth', 1.5); hold on;
plot(udds_time, Vt_est, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('V_t (V)');
title('Measured vs Estimated Terminal Voltage during UDDS Cycle');
legend('Measured V_t', 'Estimated V_t');
grid on;

% SOC 추정 함수
function [SOC_est, V1_est, Vt_est, P] = soc_estimation(SOC_est, V1_est, Vt_true, ik, Config, P, soc_values, ocv_values)
    Q = [1e-5 0; 0 1e-5]; % Process noise covariance
    R = 1e-2; % Measurement noise covariance
    
    % Prediction step
    SOC_pred = SOC_est - (Config.dt / Config.cap) * Config.coulomb_efficient * ik;
    V1_pred = exp(-Config.dt / (Config.R1 * Config.C1)) * V1_est + (1 - exp(-Config.dt / (Config.R1 * Config.C1))) * ik * Config.R1;
    X_pred = [SOC_pred; V1_pred];
    
    % State transition matrix
    A = [1, 0; 0, exp(-Config.dt / (Config.R1 * Config.C1))];
    
    % Process covariance update
    P = A * P * A' + Q;
    
    % Measurement prediction
    Vt_pred = interp1(soc_values, ocv_values, SOC_pred, 'linear', 'extrap') - Config.R1 * V1_pred - Config.R0 * ik;
    
    % Measurement update
    K = P * [1; 0] / (P(1, 1) + R); % Kalman gain
    z = Vt_true - Vt_pred; % Measurement residual
    X_pred = X_pred + K * z;
    P = (eye(2) - K * [1, 0]) * P;
    
    % Save estimates
    SOC_est = X_pred(1);
    V1_est = X_pred(2);
    Vt_est = Vt_pred + K(1) * z; % Estimated Vt
end

