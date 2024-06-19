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

% Initialize SOC estimation
SOC_est = zeros(length(udds_current), 1);
SOC_est(1) = 1; % 초기 SOC는 100%로 가정

% EKF 초기화
P = eye(2); % Initial estimation error covariance

% % Preallocate arrays for V1 and Vt
% V1_est = zeros(length(udds_current), 1);
% Vt_est = zeros(length(udds_current), 1);
% V1_est(1) = 0;
% Vt_est(1) = udds_voltage(1);

% Initialize intVar structure
intVar.SOC = zeros(length(udds_current), 1);
intVar.V1 = zeros(length(udds_current), 1);
intVar.Vt = zeros(length(udds_current), 1);

% Set initial values
intVar.SOC(1) = 1;
intVar.V1(1) = 0;
intVar.Vt(1) = udds_voltage(1);

% Configuration parameters
Config.dt = mean(diff(udds_time)); % 평균 시간 간격
Config.ik = udds_current(1); % 초기 전류
Config.R0 = 0.001884314;
Config.R1 = 0.045801322;
Config.C1 = 4846.080679;
Config.cap = (2.9950 * 3600) / 100; % nominal capacity [Ah]
Config.coulomb_efficient = 1;

% Simulation loop
for k = 2:length(udds_current)
    Config.ik = udds_current(k);
    
    % True SOC and Vt calculation
    [intVar.SOC(k), intVar.Vt(k)] = battery_model(intVar.SOC(k-1), intVar.V1(k-1), Config.ik, Config);
    
    % SOC estimation using EKF
    [intVar.SOC(k), intVar.V1(k), intVar.Vt(k), P] = soc_estimation(intVar.SOC(k-1), intVar.V1(k-1), udds_voltage(k), Config.ik, Config, P, soc_values, ocv_values);
end

% Plot SOC
figure;
plot(udds_time, intVar.SOC, 'b', 'LineWidth', 1.5); hold on;
xlabel('Time (s)');
ylabel('SOC');
title('Estimated SOC during UDDS Cycle');
grid on;

% Plot Vt
figure;
plot(udds_time, udds_voltage, 'b', 'LineWidth', 1.5); hold on;
plot(udds_time, intVar.Vt, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('V_t (V)');
title('Measured vs Estimated Terminal Voltage during UDDS Cycle');
legend('Measured V_t', 'Estimated V_t');
grid on;

% OCV 함수
function v_ocv = ocv_soc(soc, soc_values, ocv_values)
    v_ocv = interp1(soc_values, ocv_values, soc, 'linear', 'extrap');
end

% 배터리 모델 함수
function [SOC_true, Vt_true] = battery_model(SOC_prev, V1_prev, ik, Config)
    load('soc_ocv.mat', 'soc_ocv');
    soc_values = soc_ocv(:, 1);
    ocv_values = soc_ocv(:, 2);
    SOC_true = SOC_prev - (Config.dt / Config.cap) * Config.coulomb_efficient * ik;
    Vt_true = ocv_soc(SOC_true, soc_values, ocv_values) - Config.R1 * V1_prev - Config.R0 * ik;
end

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


