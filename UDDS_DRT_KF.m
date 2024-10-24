clc; clear; close all;

%% 1. 저장된 데이터 로드
% gamma_data.mat 및 soc_ocv_data.mat 파일 로드
load('gamma_data.mat', 'gamma_sorted', 'soc_sorted', 'theta_discrete', 'R0_est_all', 'soc_mid_all');
load('soc_ocv_data.mat', 'soc_values', 'ocv_values');

% tau_discrete 재정의
tau_discrete = exp(theta_discrete);

%% 2. 주행 데이터 로드
load('udds_data.mat');  % 'udds_data' 구조체를 로드합니다.

%% 3. 칼만 필터를 적용할 사이클 선택
cycle_idx = 1;  % 원하는 사이클 번호로 변경 가능
if cycle_idx > length(udds_data)
    error('선택한 사이클 번호가 데이터 범위를 벗어났습니다.');
end

% 선택한 사이클의 데이터 추출
time = udds_data(cycle_idx).t;
current_measured = udds_data(cycle_idx).I;  % 측정된 전류 (양수: 방전, 음수: 충전)
voltage_measured = udds_data(cycle_idx).V;  % 측정된 전압
SOC_true = udds_data(cycle_idx).SOC;  % 실제 SOC (있는 경우)

%% 4. 칼만 필터 설정
% 시간 간격 계산
dt = [0; diff(time)];
if dt(1) == 0
    dt(1) = dt(2);
end

% RC 소자의 개수
num_RC = length(theta_discrete);

% 초기화
num_samples = length(time);
voltage_estimated = zeros(num_samples, 1);
SOC_estimated = zeros(num_samples, 1);  % SOC 추정값 초기화

% 초기 SOC 설정 (실제 SOC가 없을 경우 0.9901로 설정)
if ~isempty(SOC_true)
    SOC_estimated(1) = SOC_true(1);
else
    SOC_estimated(1) = 0.9901;  % 초기 SOC를 99.01%로 가정
end

% 칼만 필터 파라미터 설정
state_dimension = 1 + num_RC;  % 상태 벡터: [SOC; V_RC_1; V_RC_2; ... ; V_RC_n]

% 잡음 공분산 행렬 설정
process_noise_cov = eye(state_dimension) * 1e-7; % 프로세스 잡음 공분산 (필요에 따라 조정)
measurement_noise_cov = 1e-16;  % 측정 잡음 공분산 (필요에 따라 조정)

% 초기 상태 추정
state_estimate = zeros(state_dimension, 1);
state_estimate(1) = SOC_estimated(1);
state_estimate(2:end) = 0;  % 초기 RC 전압은 0으로 가정

% 초기 공분산 행렬
state_covariance = eye(state_dimension) * 1e-3;

%% 5. 칼만 필터 적용
for k = 1:num_samples-1
    % 현재 SOC에 대한 gamma 값 보간
    gamma_current = interp1(soc_sorted, gamma_sorted, SOC_estimated(k), 'linear', 'extrap');  % (1 x num_RC)
    
    % 각 RC 소자에 대한 R_i 계산
    delta_theta = theta_discrete(2) - theta_discrete(1);  % 스칼라
    R_i = gamma_current * delta_theta;  % (1 x num_RC)
    
    % C_i 계산
    C_i = tau_discrete ./ R_i;  % (1 x num_RC)
    
    % 상태 전이 행렬 A 및 제어 입력 행렬 B 계산
    A = eye(state_dimension);
    for i = 1:num_RC
        A(i+1, i+1) = exp(-dt(k)/tau_discrete(i));  % RC 전압의 상태 전이
    end
    
    % 제어 입력 행렬 B 계산
    B = zeros(state_dimension, 1);
    B(1) = dt(k)/3600; 
    for i = 1:num_RC
        B(i+1) = R_i(i) * (1 - exp(-dt(k)/tau_discrete(i)));
    end
    
    % 시스템의 제어 입력 (전류)
    control_input = current_measured(k);
    
    % 예측 단계
    state_prediction = A * state_estimate + B * control_input;
    state_covariance_prediction = A * state_covariance * A' + process_noise_cov;
    
    % OCV 계산
    OCV_predicted = interp1(soc_values, ocv_values, state_prediction(1), 'linear', 'extrap');
    
    % 관측 행렬 H 계산
    delta_SOC = 1e-5;
    OCV_plus = interp1(soc_values, ocv_values, state_prediction(1) + delta_SOC, 'linear', 'extrap');
    OCV_minus = interp1(soc_values, ocv_values, state_prediction(1) - delta_SOC, 'linear', 'extrap');
    dOCV_dSOC = (OCV_plus - OCV_minus) / (2 * delta_SOC);
    
    H = zeros(1, state_dimension);
    H(1) = dOCV_dSOC;  % SOC에 대한 OCV의 미분
    H(2:end) = 1;  % 각 RC 전압의 계수
    
    % 예측된 전압
    V_predicted = OCV_predicted + R0_est_all(cycle_idx) * control_input + sum(state_prediction(2:end));
    voltage_estimated(k) = V_predicted;
    
    % 칼만 이득 계산
    S = H * state_covariance_prediction * H' + measurement_noise_cov;  % 스케일링 팩터
    kalman_gain = (state_covariance_prediction * H') / S;  % (state_dimension x 1)
    
    % 측정 업데이트
    measurement_residual = voltage_measured(k) - V_predicted;  % 잔차
    state_estimate = state_prediction + kalman_gain * measurement_residual;
    state_covariance = (eye(state_dimension) - kalman_gain * H) * state_covariance_prediction;
    
    % 상태 저장
    SOC_estimated(k+1) = state_estimate(1);
end

% 마지막 전압 추정
V_predicted_final = interp1(soc_values, ocv_values, SOC_estimated(end), 'linear', 'extrap') + ...
                    R0_est_all(cycle_idx) * current_measured(end) + sum(state_estimate(2:end));
voltage_estimated(end) = V_predicted_final;

%% 6. 결과 시각화
figure;
subplot(3,1,1);
plot(time, SOC_true, 'b', 'LineWidth', 1.5); hold on;
plot(time, SOC_estimated, 'r--', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('SOC');
title('실제 SOC vs 추정 SOC');
legend('실제 SOC', '추정 SOC');
grid on;

subplot(3,1,2);
plot(time, voltage_measured, 'b', 'LineWidth', 1.5); hold on;
plot(time, voltage_estimated, 'r--', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('전압 (V)');
title('측정 전압 vs 추정 전압');
legend('측정 전압', '추정 전압');
grid on;

subplot(3,1,3);
plot(time, current_measured, 'k', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('전류 (A)');
title('전류 프로파일');
grid on;

%% 7. 추가 시각화 및 결과 저장 (필요 시)
% SOC 추정 오류 계산
if ~isempty(SOC_true)
    SOC_error = SOC_true - SOC_estimated;
    SOC_RMSE = sqrt(mean(SOC_error.^2));
    fprintf('SOC RMSE: %.4f%%\n', SOC_RMSE * 100);
    
    % SOC 추정 오류 시각화
    figure;
    plot(time, SOC_error, 'g', 'LineWidth', 1.5);
    xlabel('시간 (s)');
    ylabel('SOC 오류');
    title('SOC 추정 오류 시간에 따른 변화');
    legend(sprintf('SOC RMSE: %.4f%%', SOC_RMSE * 100));
    grid on;
end

% 결과 행렬 저장
result_matrix = [SOC_true, SOC_estimated, state_estimate(2:end)'];
column_names = {'True_SOC', 'Estimated_SOC', 'V_RC_1', 'V_RC_2', 'V_RC_3', ... % 필요한 만큼 추가
                'V_RC_n'};
%save('kalman_filter_result.mat', 'result_matrix', 'column_names');

% 전압 잔차 계산 및 시각화
voltage_residuals = voltage_measured - voltage_estimated;
R_calculated = var(voltage_residuals);
fprintf('Calculated sensor noise covariance R: %.4f\n', R_calculated);

figure;
plot(time, voltage_residuals, 'k', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('잔차 (V)');
title('측정 전압과 추정 전압의 잔차');
grid on;

