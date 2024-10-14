clear; clc; close all;

%% 데이터 로드
load('C:\Users\USER\Desktop\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
udds_current = meas.Current;  % 전류 데이터
udds_voltage = meas.Voltage;  % 전압 데이터
udds_time = meas.Time;        % 시간 데이터

% SOC-OCV 데이터 로드
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);  % SOC 값
ocv_values = soc_ocv(:, 2);  % OCV 값

%% 초기 설정
initial_soc = 0.9901;          % 초기 SOC 값
Config.dt = mean(diff(udds_time));  % 시간 간격
Config.cap = 2.90;                  % 배터리 용량 (Ah)
Config.coulomb_efficient = 1;       % 쿨롱 효율

%% DRT 설정
n = 21;                         % RC 회로의 수
lambda = 0.0409;                % 정규화 파라미터
mu = 10;                        % 저항 분포의 평균
sigma = 5;                      % 저항 분포의 표준 편차
tau_discrete = linspace(udds_time(1), 2500, n);  % 이산 타임 상수

% 초기 저항 분포 (정규 분포 기반)
R_discrete_true = normpdf(tau_discrete, mu, sigma);
R_discrete_true = R_discrete_true / max(R_discrete_true);  % 최대값 1로 정규화

% 초기 OCV 추정값
OCV_est = interp1(soc_values, ocv_values, initial_soc, 'linear', 'extrap');

%% 전압 및 전류 초기화
V_est = zeros(1, length(udds_time));  % n-RC 모델을 통한 전압 계산
R0 = 0.025426;                            % 초기 저항값
V_RC = zeros(n, length(udds_time));   % 각 RC 소자의 전압 저장 공간

%% 전압 계산 (첫 번째 시간 스텝)
for i = 1:n
    V_RC(i, 1) = udds_current(1) * R_discrete_true(i) * (1 - exp(-Config.dt / tau_discrete(i)));
end
V_est(1) = OCV_est + R0 * udds_current(1) + sum(V_RC(:, 1));

%% 전압 계산 (t > 1)
for k = 2:length(udds_time)
    for i = 1:n
        V_RC(i, k) = exp(-Config.dt / tau_discrete(i)) * V_RC(i, k-1) + ...
                     R_discrete_true(i) * (1 - exp(-Config.dt / tau_discrete(i))) * udds_current(k);
    end
    V_est(k) = OCV_est + R0 * udds_current(k) + sum(V_RC(:, k));
end

%% 정규화 행렬 L 생성 (1차 차분)
L = zeros(n-1, n);
for i = 1:n-1
    L(i, i) = -1;
    L(i, i+1) = 1;
end

%% W 행렬 구성
W = zeros(length(udds_time), n);  % W 행렬 초기화
for k = 1:length(udds_time)
    for i = 1:n
        if k == 1
            W(k, i) = udds_current(k) * (1 - exp(-Config.dt / tau_discrete(i)));
        else
            W(k, i) = exp(-Config.dt / tau_discrete(i)) * W(k-1, i) + ...
                      udds_current(k) * (1 - exp(-Config.dt / tau_discrete(i)));
        end
    end
end

%% 정규화된 DRT 추정 (정규화 적용)
R_estimated = (W' * W + lambda * (L' * L)) \ (W' * udds_voltage); 
R_estimated(R_estimated < 0) = 0;  % 음수 값은 0으로 설정

%% 결과 출력
figure;
%plot(tau_discrete, R_discrete_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True DRT');
%hold on;
plot(tau_discrete, R_estimated, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Estimated DRT');
xlabel('\tau (Time Constant)');
ylabel('R (Resistance)');
legend('Location', 'Best');
title('True DRT vs Estimated DRT with UDDS Data');
grid on;