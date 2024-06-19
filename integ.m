clc; clear; close all;

% Load UDDS data
load('C:\Users\USER\Desktop\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
udds_current = meas.Current; % UDDS 전류 데이터
udds_voltage = meas.Voltage; % UDDS 전압 데이터
udds_time = meas.Time; % UDDS 시간 데이터

% Load SOC-OCV structure
load('soc_ocv.mat', 'soc_ocv'); %c/20에서 가져온 soc-ocv lookup table
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

% Remove duplicate OCV values
[unique_ocv_values, unique_idx] = unique(ocv_values);
unique_soc_values = soc_values(unique_idx);

% Interpolate SOC for the first UDDS voltage value
initial_voltage = udds_voltage(1);
initial_soc = interp1(unique_ocv_values, unique_soc_values, initial_voltage, 'linear', 'extrap');

fprintf('Initial voltage: %.2f V corresponds to SOC: %.2f%%\n', initial_voltage, initial_soc * 100);

% Initialize SOC estimation
SOC_est = zeros(length(udds_current), 1);
SOC_est(1) = initial_soc; % 초기 SOC는 100%로 가정 % 완충상태 가정 (udds voltage 4.18v) 

% Initialize true SOC using coulomb counting
true_SOC = zeros(length(udds_current), 1);
true_SOC(1) = initial_soc; % 초기 SOC는 100%로 가정

% EKF 초기화
P = eye(2); % Initial estimation error covariance

% Preallocate arrays for V1 and Vt
V1_est = zeros(length(udds_current), 1);
Vt_est = zeros(length(udds_current), 1);
V1_est(1) = udds_current(1) * Config.R1 * (1 - exp(-Config.dt / (Config.R1 * Config.C1))); % 초기 V1 값 (전 상태의 v1값이 없기때문에 앞의 항이 zero)
Vt_est(1) = udds_voltage(1); % 초기 측정 전압 = estimated 전압 

% Simulation loop
for k = 2:length(udds_current)
    Config.ik = udds_current(k); % 샘플링 시간동안 흐르는 전류 

    % True SOC calculation using coulomb counting
    delta_t = udds_time(k) - udds_time(k-1); % dt 써도 되는데, 정확한 값을 위하여... (비교해보니까 거의 비슷하긴 함) 
    true_SOC(k) = true_SOC(k-1) + (udds_current(k) * delta_t) / (Config.cap * 3600); % 실제 soc = 전 soc + i * dt/q_total (sampling 시간 동안 흐르는 전류)

    % SOC estimation using EKF
    [SOC_est(k), V1_est(k), Vt_est(k), P] = soc_estimation(SOC_est(k-1), V1_est(k-1), udds_voltage(k), udds_current(k), Config, P, soc_values, ocv_values); % 추정 코드
end

% Plot SOCrm
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
    
    % Prediction step (상태방정식)
    SOC_pred = SOC_est + (Config.dt / Config.cap) * Config.coulomb_efficient * ik;
    V1_pred = exp(-Config.dt / (Config.R1 * Config.C1)) * V1_est + (1 - exp(-Config.dt / (Config.R1 * Config.C1))) * ik * Config.R1;
    X_pred = [SOC_pred; V1_pred];
    
    % State transition matrix
    A = [1, 0;
         0, exp(-Config.dt / (Config.R1 * Config.C1))];
    
    % Process covariance update
    P = A * P * A' + Q;
    
    % Measurement prediction (OCV = 전류적산법으로 얻은 SOC로부터 lookup table 통과시켜 얻음)
    Vt_pred = interp1(soc_values, ocv_values, SOC_pred, 'linear', 'extrap') - Config.R1 * V1_pred - Config.R0 * ik;
    
    % H 행렬 계산
    dOCV_dSOC = (interp1(soc_values, ocv_values, SOC_pred + 1e-5, 'linear', 'extrap') - ...
                 interp1(soc_values, ocv_values, SOC_pred, 'linear', 'extrap')) / 1e-5;
    H = [dOCV_dSOC, -1];
    
    % Measurement update
    K = P * H' / (H * P * H' + R); % Kalman gain
    z = Vt_true - Vt_pred; % Measurement residual
    X_pred = X_pred + K * z; % 상태 변수 업데이트
    P = (eye(2) - K * H) * P; % 공분산 업데이트
    
    % Save estimates
    SOC_est = X_pred(1);
    V1_est = X_pred(2);
    Vt_est = Vt_pred + K(1) * z; % Estimated Vt 업데이트
end

