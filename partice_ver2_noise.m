clc; clear; close all;

%% 1. 데이터 로드

% 데이터 파일 경로를 사용자 환경에 맞게 수정하세요.
load('G:\공유 드라이브\BSL_Data3\Driving cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
load('C:\Users\deu04\OneDrive\바탕 화면\ECM_github\kalman_new\Real Data\udds_data.mat');

% 개별 트립 데이터
trip1_I = udds_data(1).I;
trip1_V = udds_data(1).V;
trip1_t = udds_data(1).t;

% 전체 트립 데이터
All_trips_I = meas.Current;
All_trips_V = meas.Voltage; 
All_trips_t = meas.Time; 

SOC_initial = 0.9901;

%% 2. Markov chain noise를 전류에 추가 및 SOC 계산

n = 21;

noise_vector = linspace(min(trip1_I), max(trip1_I), n);

P = zeros(n);

% 표준 편차 설정 (필요에 따라 조정 가능)
sigma = (max(noise_vector) - min(noise_vector)) / 10;

% 각 행에 대해 전이 확률 계산
for i = 1:n
    % 현재 상태에서 가능한 모든 다음 상태에 대한 확률 계산
    probabilities = normpdf(noise_vector, noise_vector(i), sigma);
    
    % 정규화하여 확률의 합이 1이 되도록 함
    P(i, :) = probabilities / sum(probabilities);
end

%% 추가: Markov 체인을 사용하여 노이즈 추가 및 상태 기록

noisy_I = zeros(size(trip1_I));
states = zeros(size(trip1_I)); % 상태를 저장할 배열
current_state = 1; % 초기 상태 설정 (예: noise_vector의 첫 번째 값)

for t = 1:length(trip1_I)
    % 현재 상태에 해당하는 노이즈 값을 추가
    noisy_I(t) = trip1_I(t) + noise_vector(current_state);
    
    % 상태 기록
    states(t) = current_state;
    
    % 다음 상태로 전이
    current_state = find(rand <= cumsum(P(current_state, :)), 1, 'first');
end

% SOC 계산
SOC_noisy = SOC_initial + cumtrapz(trip1_t/3600, noisy_I)/2.9; % 예시: 누적 전류에 비례하여 SOC 감소
SOC_true  = SOC_initial + cumtrapz(trip1_t/3600, trip1_I)/2.9 ;

%% 결과 시각화
figure;

subplot(4,1,1);
plot(trip1_t, trip1_I);
title('Original Current');
xlabel('Time (s)');
ylabel('Current (A)');

subplot(4,1,2);
plot(trip1_t, noisy_I);
title('Noisy Current with Markov Chain');
xlabel('Time (s)');
ylabel('Current (A)');

subplot(4,1,3);
plot(trip1_t, SOC_noisy);
hold on
plot(trip1_t, SOC_true);
title('State of Charge (SOC)');
xlabel('Time (s)');
ylabel('SOC (%)');
legend('Noisy SOC', 'True SOC');

%% 3. Markov 체인 노이즈 확인을 위한 시각화

% 3.1. 상태 변화 시각화 (시간에 따른 상태 변화)
subplot(4,1,4);
plot(trip1_t, noise_vector(states));
title('Noise States Over Time');
xlabel('Time (s)');
ylabel('Noise Value (A)');

% 3.2. 상태 분포 히스토그램
figure;
histogram(noise_vector(states), 'Normalization', 'probability');
title('Histogram of Noise States');
xlabel('Noise Value (A)');
ylabel('Probability');
grid on;

% 3.3. 상태 전이 다이어그램
figure;
G = digraph(P);
figure;
plot(G, 'Layout', 'force', 'EdgeLabel', round(P,2));
title('State Transition Diagram');
xlabel('States');
ylabel('States');

% 추가: 상태 클러스터링 확인을 위한 2D 히스토그램 또는 다른 시각화 방법도 고려할 수 있습니다.

% 
