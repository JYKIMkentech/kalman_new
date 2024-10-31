clc; clear; close all;

%% 1. 데이터 로드
% udds_data.mat 파일이 현재 작업 디렉토리에 있어야 합니다.

load('G:\공유 드라이브\BSL_Data3\Driving cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');

%load('udds_data.mat');

% I = udds_data(1).I;        % 전류 (A)
% V = udds_data(1).V;        % 전압 (V)
% t = udds_data(1).t;        % 시간 (초)
% CC_SOC = udds_data(1).SOC; % Coulomb 카운팅을 통한 SOC (%)

I = meas.Current; % UDDS 전류 데이터
V = meas.Voltage; % UDDS 전압 데이터
t = meas.Time; % UDDS 시간 데이터



%% 2. 노이즈 상태 벡터 구성

% 2.1 전류 데이터의 평균 절대값 계산
I_mean = mean(abs(I)); % 평균 전류의 절대값 (A)

% 2.2 평균 전류의 10% 계산하여 노이즈 최대값 설정
n_max = I_mean * 0.1; % 노이즈 최대값 (A)

% 2.3 노이즈 상태 벡터 정의
noise_states = linspace(-n_max, n_max, 5); % 노이즈 상태 (-n_max ~ n_max 사이의 5개 값)

% 2.5 노이즈의 표준편차 설정 (노이즈 상태 내의 작은 변동)
sigma = n_max * 0.1; % 노이즈 최대값의 10%

%% 3. 전이 확률 행렬 정의

% 전이 확률 행렬 정의 (5x5 행렬)
P = [0.7, 0.15, 0.1, 0.05, 0;
     0.15, 0.7, 0.1, 0.05, 0;
     0.1, 0.1, 0.6, 0.1, 0.1;
     0.05, 0.1, 0.1, 0.7, 0.05;
     0, 0.05, 0.1, 0.15, 0.7];

% 전이 확률 행렬이 유효한지 확인 (각 행의 합이 1인지 확인)
assert(all(abs(sum(P,2) - 1) < 1e-6), '전이 확률 행렬의 각 행의 합이 1이 아닙니다.');

%% 4. 노이즈 생성

% 4.1 노이즈 벡터 초기화
n = zeros(size(I));

% 4.2 초기 상태 설정 (노이즈 상태 벡터의 인덱스)
current_state = 3; % 초기 상태는 0 A (인덱스 3)

% 4.3 노이즈 생성 루프 (PDF 기반 샘플링)
% 추가: 상태 시퀀스를 저장하기 위한 벡터 초기화
states_sequence = zeros(size(I)); 
states_sequence(1) = noise_states(current_state); % 첫 번째 상태 저장

for k = 1:length(I)
    % 현재 상태의 전이 확률 추출
    transition_prob = P(current_state, :);
    
    % 다음 상태 샘플링 (PDF 기반)
    next_state = randsample(noise_states, 1, true, transition_prob);
    
    % 노이즈 계산 (상태 값에 작은 변동 추가)
    n(k) = next_state + sigma * randn();
    
    % 다음 상태의 인덱스 찾기
    current_state = find(noise_states == next_state, 1);
    
    % 상태 시퀀스 저장
    states_sequence(k) = next_state;
end

% 4.4 노이즈가 적용된 전류 생성
I_noisy = I + n;

%% 5. SOC 계산

% 5.1 파라미터 설정
C_nom = 2.9 * 3600;          % 배터리의 정격 용량 (Ah)
t_hours = t ;           % 시간을 시간 단위로 변환 (초 -> 시간)
SOC_initial = 0.9901;         % 초기 SOC (비율, 1은 100%)

% 5.2 SOC 계산 (cumtrapz 함수 사용)
% 충전 시 SOC가 증가하고 방전 시 SOC가 감소하도록 수식을 수정
% 방전 시 전류가 양수일 경우 SOC 감소, 충전 시 음수일 경우 SOC 증가
SOC_original = SOC_initial + cumtrapz(t_hours, I) / C_nom;
SOC_noisy = SOC_initial + cumtrapz(t_hours, I_noisy) / C_nom;

%% 6. 결과 시각화

% 6.1 전류 비교
figure;
subplot(3,2,1); % 첫 번째 subplot을 첫 번째 위치에 설정
plot(t, I, 'b', 'LineWidth', 1.5); hold on;
plot(t, I_noisy, 'r', 'LineWidth', 1.5);
xlabel('시간 (초)');
ylabel('전류 (A)');
title('원래 전류 vs. Markov noise가 적용된 전류');
legend('원래 전류', 'Markov noise 적용 전류');
grid on;

% 6.2 SOC 비교
subplot(3,2,2); % 두 번째 subplot을 두 번째 위치에 설정
hold on
plot(t, SOC_original , 'b-', 'LineWidth', 1.5); % 비율을 퍼센트로 변환
plot(t, SOC_noisy , 'r--', 'LineWidth', 1.5);
xlabel('시간 (초)');
ylabel('SOC (%)');
title('True SOC vs. Markov noise가 적용된 CC SOC');
legend('True SOC', 'Markov noise 적용한 CC SOC');
grid on;

% 6.3 SOC 차이 그래프 추가
subplot(3,2,3); % 세 번째 subplot을 세 번째 위치에 설정
SOC_difference = (SOC_original - SOC_noisy ) * 100; % 퍼센트 단위로 변환
plot(t, SOC_difference, 'm', 'LineWidth', 1.5);
xlabel('시간 (초)');
ylabel('SOC 차이 (%)');
title('Markov noise가 적용된 SOC와 원래 SOC의 차이');
legend('SOC 차이');
grid on;

%% 7. Markov 노이즈 검증

% 7.1 노이즈 상태 시퀀스 시각화
subplot(3,2,4);
plot(t, states_sequence, 'LineWidth', 1.5);
xlabel('시간 (초)');
ylabel('노이즈 상태 (A)');
title('Markov 노이즈 상태 시퀀스');
grid on;

% 7.2 상태 전이 히스토그램 및 전이 확률 비교
% 상태 전이 횟수 계산
transition_counts = zeros(length(noise_states));
for k = 2:length(states_sequence)
    from = find(noise_states == states_sequence(k-1));
    to = find(noise_states == states_sequence(k));
    transition_counts(from, to) = transition_counts(from, to) + 1;
end

% 전이 확률 계산
transition_prob_estimated = transition_counts ./ sum(transition_counts, 2);

% 전이 확률 비교를 위한 subplot
subplot(3,2,5);
imagesc(transition_prob_estimated); % 전이 확률 행렬 시각화
colorbar;
xlabel('다음 상태');
ylabel('현재 상태');
title('추정된 전이 확률 행렬');
set(gca, 'XTick', 1:length(noise_states), 'XTickLabel', noise_states);
set(gca, 'YTick', 1:length(noise_states), 'YTickLabel', noise_states);

% 정의된 전이 확률 행렬 P와 추정된 전이 확률 행렬 비교
disp('정의된 전이 확률 행렬 P:');
disp(P);
disp('추정된 전이 확률 행렬:');
disp(transition_prob_estimated);

% 7.3 상태 빈도수 히스토그램
subplot(3,2,6);
histogram(states_sequence, noise_states, 'Normalization', 'probability');
xlabel('노이즈 상태 (A)');
ylabel('빈도');
title('노이즈 상태 빈도수 히스토그램');
grid on;

% 전체 그림 제목 설정
sgtitle('Markov 노이즈 적용 및 검증 결과');

%% 8. 추가 설명
% 위의 플롯들을 통해 다음을 확인할 수 있습니다:
% - 원래 전류와 노이즈가 적용된 전류의 비교를 통해 노이즈가 올바르게 추가되었는지 확인.
% - SOC 비교 및 SOC 차이 그래프를 통해 노이즈가 SOC 계산에 미치는 영향을 평가.
% - 노이즈 상태 시퀀스 시각화를 통해 상태들이 Markov 체인에 따라 클러스터링되어 있는지 확인.
% - 전이 확률 행렬을 비교하여 생성된 노이즈의 전이 확률이 정의한 P와 일치하는지 검증.
% - 상태 빈도수 히스토그램을 통해 각 상태가 예상대로 나타나는지 확인.


