clear; clc; close all;

%% 1. 저장된 데이터 로드
% gamma_data.mat 및 soc_ocv_data.mat 파일 로드
load('gamma_data.mat', 'gamma_sorted', 'soc_sorted', 'theta_discrete', 'R0_est_all', 'soc_mid_all');
load('soc_ocv_data.mat', 'soc_values', 'ocv_values');

% tau_discrete 재정의
tau_discrete = exp(theta_discrete);  % tau_discrete 재정의

%% 2. 주행 데이터 로드
load('udds_data.mat');  % 'udds_data' 구조체를 로드합니다.

%% 3. 칼만 필터를 적용할 사이클 선택
cycle_idx = 1;  % 원하는 사이클 번호로 변경 가능
if cycle_idx > length(udds_data)
    error('선택한 사이클 번호가 데이터 범위를 벗어났습니다.');
end

% 선택한 사이클의 데이터 추출
t = udds_data(cycle_idx).t;
I_meas = udds_data(cycle_idx).I;  % 측정된 전류 (양수: 방전, 음수: 충전)
V_meas = udds_data(cycle_idx).V;  % 측정된 전압
SOC_true = udds_data(cycle_idx).SOC;  % 실제 SOC (있는 경우)

%% 4. 칼만 필터 설정
% 시간 간격 계산
dt = [0; diff(t)];
if dt(1) == 0
    dt(1) = dt(2);
end

% RC 소자의 개수
n_RC = length(theta_discrete);

% 초기화
N = length(t);
V_est = zeros(N, 1);
SOC_est = zeros(N, 1);  % SOC 추정값 초기화

% 초기 SOC 설정 (실제 SOC가 없을 경우 0.5 등으로 설정 가능)
if ~isempty(SOC_true)
    SOC_est(1) = SOC_true(1);
else
    SOC_est(1) = 0.5;  % 초기 SOC를 50%로 가정
end

% 칼만 필터 파라미터 설정
% 상태 벡터: [SOC; V_RC_1; V_RC_2; ... ; V_RC_n]
state_dim = 1 + n_RC;

% 상태 전이 행렬 A와 제어 입력 행렬 B는 시간에 따라 변하므로 루프 내에서 설정
% 관측 행렬 H는 전압 측정을 위한 행렬로 설정

% 잡음 공분산 행렬 설정
Q = eye(state_dim) * 1e-5;  % 프로세스 잡음 공분산 (필요에 따라 조정)
R_kf = 1e-3;  % 측정 잡음 공분산 (필요에 따라 조정)

% 초기 상태 추정
x_est = zeros(state_dim, 1);
x_est(1) = SOC_est(1);
% 초기 RC 전압은 0으로 가정
x_est(2:end) = 0;

% 초기 공분산 행렬
P = eye(state_dim) * 1e-3;

%% 5. 칼만 필터 적용
for k = 1:N-1
    % 현재 SOC에 대한 gamma 값 보간
    gamma_k = interp1(soc_sorted, gamma_sorted, SOC_est(k), 'linear', 'extrap');  % (1 x n_RC)
    
    % 각 RC 소자에 대한 R_i 계산
    delta_theta = theta_discrete(2) - theta_discrete(1);  % 스칼라
    R_i = gamma_k * delta_theta;  % (1 x n_RC)
    
    % C_i 계산
    C_i = tau_discrete ./ R_i;  % (1 x n_RC)
    
    % 상태 전이 행렬 A 및 제어 입력 행렬 B 계산
    % A는 (state_dim x state_dim) 행렬
    A = eye(state_dim);
    for i = 1:n_RC
        A(1, 1) = 1;  % SOC는 상태 유지
        A(i+1, i+1) = exp(-dt(k)/tau_discrete(i));  % RC 전압의 상태 전이
    end
    
    % 제어 입력 행렬 B 계산
    B = zeros(state_dim, 1);
    B(1) = -dt(k)/3600;  % SOC 변화: C / Q = I * dt / Q (Q를 1으로 가정)
    for i = 1:n_RC
        B(i+1) = R_i(i) * (1 - exp(-dt(k)/tau_discrete(i)));
    end
    
    % 시스템의 제어 입력 (전류)
    u = I_meas(k);
    
    % 예측 단계
    x_pred = A * x_est + B * u;
    P_pred = A * P * A' + Q;
    
    % OCV 계산
    OCV_pred = interp1(soc_values, ocv_values, x_pred(1), 'linear', 'extrap');
    
    % 관측 행렬 H 계산
    % H는 (1 x state_dim) 행렬
    % V_est = OCV + R0 * I + sum(V_RC)
    % dOCV/dSOC 계산 (수치 미분)
    delta_SOC = 1e-5;
    OCV_plus = interp1(soc_values, ocv_values, x_pred(1) + delta_SOC, 'linear', 'extrap');
    OCV_minus = interp1(soc_values, ocv_values, x_pred(1) - delta_SOC, 'linear', 'extrap');
    dOCV_dSOC = (OCV_plus - OCV_minus) / (2 * delta_SOC);
    
    H = zeros(1, state_dim);
    H(1) = dOCV_dSOC;  % SOC에 대한 OCV의 미분
    H(2:end) = 1;  % 각 RC 전압의 계수
    
    % 예측된 전압
    V_pred = OCV_pred + R0_est_all(cycle_idx) * u + sum(x_pred(2:end));
    V_est(k) = V_pred;
    
    % 칼만 이득 계산
    S = H * P_pred * H' + R_kf;  % 스케일링 팩터
    K = (P_pred * H') / S;  % (state_dim x 1)
    
    % 측정 업데이트
    y = V_meas(k) - V_pred;  % 잔차
    x_est = x_pred + K * y;
    P = (eye(state_dim) - K * H) * P_pred;
    
    % 상태 저장
    SOC_est(k+1) = x_est(1);
end

% 마지막 전압 추정
V_est(N) = interp1(soc_values, ocv_values, SOC_est(N), 'linear', 'extrap') + R0_est_all(cycle_idx) * I_meas(N) + sum(x_est(2:end));

%% 6. 결과 시각화
figure;
subplot(3,1,1);
plot(t, SOC_true, 'b', 'LineWidth', 1.5); hold on;
plot(t, SOC_est, 'r--', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('SOC');
title('실제 SOC vs 추정 SOC');
legend('실제 SOC', '추정 SOC');
grid on;

subplot(3,1,2);
plot(t, V_meas, 'b', 'LineWidth', 1.5); hold on;
plot(t, V_est, 'r--', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('전압 (V)');
title('측정 전압 vs 추정 전압');
legend('측정 전압', '추정 전압');
grid on;

subplot(3,1,3);
plot(t, I_meas, 'k', 'LineWidth', 1.5);
xlabel('시간 (s)');
ylabel('전류 (A)');
title('전류 프로파일');
grid on;

