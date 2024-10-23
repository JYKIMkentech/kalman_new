clear; clc; close all;

%% 1. UDDS 주행 데이터 로드
% UDDS 주행 데이터를 로드합니다.
load('C:\Users\deu04\OneDrive\바탕 화면\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
load('udds_data.mat');
udds_current = meas.Current;  % 전류 데이터 (A)
udds_voltage = meas.Voltage;  % 전압 데이터 (V)
udds_time = meas.Time;        % 시간 데이터 (s)

% 시간 벡터가 0에서 시작하고 연속적이도록 보정합니다.
udds_time = udds_time - udds_time(1);

%% 2. SOC-OCV 데이터 로드
% SOC-OCV 데이터를 로드합니다.
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);  % SOC 값
ocv_values = soc_ocv(:, 2);  % OCV 값

%% 3. SOC 계산
% 초기 SOC 설정
SOC_initial = 0.9901;  % 초기 SOC를 100%로 가정합니다.

% 전류 데이터의 시간 간격(delta_t) 계산
delta_t = [0; diff(udds_time)];  % 각 측정 사이의 시간 간격

% 배터리 용량 (Ah를 Coulomb로 변환)
Q_battery = 2.9 * 3600;  % 배터리 용량 (2.9Ah)

% 시간에 따른 SOC 계산
udds_SOC = SOC_initial + cumsum(udds_current .* delta_t) / Q_battery;

%% Plot

plot(udds_time,udds_current);
hold on
plot(udds_time,udds_SOC);
xlabel('time')
ylabel('current')
yyaxis right
ylabel('soc')


