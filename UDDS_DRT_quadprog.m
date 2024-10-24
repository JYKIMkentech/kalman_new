clear; clc; close all;

%% 1. UDDS 주행 데이터 로드
% UDDS 주행 데이터를 로드합니다.
load('udds_data.mat');  % 'udds_data' 구조체를 로드합니다.

%% 2. SOC-OCV 데이터 로드
% SOC-OCV 데이터를 로드합니다.
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);  % SOC 값
ocv_values = soc_ocv(:, 2);  % OCV 값

%% 3. DRT 추정에 필요한 파라미터 설정
n = 40;  % 이산 요소의 개수
tau_min = 0.1;     % 최소 시간 상수 (초)
tau_max = 1100;    % 최대 시간 상수 (초)

% Theta 및 tau 값 계산
theta_min = log(tau_min);
theta_max = log(tau_max);
theta_discrete = linspace(theta_min, theta_max, n);
tau_discrete = exp(theta_discrete);

% Delta theta 계산
delta_theta = theta_discrete(2) - theta_discrete(1);

% 정규화 파라미터
lambda = 0.153;  % 필요에 따라 조정 가능

% Gamma에 대한 1차 차분 행렬 L_gamma 생성
L_gamma = zeros(n-1, n);
for i = 1:n-1
    L_gamma(i, i) = -1;
    L_gamma(i, i+1) = 1;
end

% R0에 대한 정규화를 피하기 위해 L_aug 생성
L_aug = [L_gamma, zeros(n-1, 1)];

%% 4. 각 트립에 대한 DRT 추정 (quadprog 사용)
num_trips = length(udds_data);

% 결과 저장을 위한 배열 사전 할당
gamma_est_all = zeros(num_trips-1, n);  % 마지막 트립 제외
R0_est_all = zeros(num_trips-1, 1);
soc_mid_all = zeros(num_trips-1, 1);  % 각 트립의 중간 SOC 저장

for s = 1:num_trips-1  % 마지막 트립은 데이터가 짧으므로 제외
    fprintf('Processing Trip %d/%d...\n', s, num_trips-1);
    
    % 현재 트립의 데이터 추출
    t = udds_data(s).t;
    ik = udds_data(s).I;
    V_sd = udds_data(s).V;
    SOC = udds_data(s).SOC;
    
    % 각 트립의 중간 SOC 계산 (중간 시간에 해당하는 SOC)
    t_mid = t(end) / 2;  % 트립의 중간 시간
    
    % t와 SOC를 고유한 t에 대해 정렬 및 중복 제거
    [t_unique, idx_unique] = unique(t);
    SOC_unique = SOC(idx_unique);
    
    soc_mid_all(s) = interp1(t_unique, SOC_unique, t_mid, 'linear', 'extrap');  % 중간 시간에 해당하는 SOC
    
    % 시간 간격 dt 계산
    delta_t = [0; diff(t)];
    dt = delta_t;
    if dt(1) == 0  % 첫 번째 dt 값이 0이면
        dt(1) = dt(2);  % 두 번째 dt 값으로 대체
    end
    
    % OCV 계산 (SOC-OCV 테이블 사용)
    ocv_over_time = interp1(soc_values, ocv_values, SOC, 'linear', 'extrap');
    
    %% 4.1 W 행렬 생성
    % W 행렬 초기화
    W = zeros(length(t), n);
    for k_idx = 1:length(t)
        for i = 1:n
            if k_idx == 1
                W(k_idx, i) = ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i))) * delta_theta;
            else
                W(k_idx, i) = W(k_idx-1, i) * exp(-dt(k_idx) / tau_discrete(i)) + ...
                              ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i))) * delta_theta;
            end
        end
    end
    
    % R0 추정을 위한 W 행렬 확장
    W_aug = [W, ik(:)];  % ik(:)는 ik를 열 벡터로 변환
    
    %% 4.2 y 벡터 생성
    % y = V_sd - OCV
    y = V_sd - ocv_over_time;
    y = y(:);  % y를 열 벡터로 변환
    
    %% 4.3 quadprog를 사용한 제약 조건 하의 추정
    % 비용 함수: 0.5 * Theta' * H * Theta + f' * Theta
    H = (W_aug' * W_aug + lambda * (L_aug' * L_aug));
    f = -W_aug' * y;
    
    % 제약 조건: Theta >= 0 (gamma와 R0는 0 이상)
    A = -eye(n+1);
    b = zeros(n+1, 1);
    
    % quadprog 옵션 설정
    options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
    
    % quadprog 실행
    [Theta_est, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], options);
    
    if exitflag ~= 1
        warning('Optimization did not converge for trip %d.', s);
    end
    
    % gamma와 R0 추정값 추출
    gamma_est = Theta_est(1:n);
    R0_est = Theta_est(n+1);
    
    % 추정값 저장
    gamma_est_all(s, :) = gamma_est';
    R0_est_all(s) = R0_est;
    
    %% 4.4 V_est 계산 (검증용)
    % V_RC 및 V_est 초기화
    V_RC = zeros(n, length(t));  % 각 요소의 전압
    V_est = zeros(length(t), 1);
    for k_idx = 1:length(t)
        if k_idx == 1
            for i = 1:n
                V_RC(i, k_idx) = gamma_est(i) * delta_theta * ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i)));
            end
        else
            for i = 1:n
                V_RC(i, k_idx) = V_RC(i, k_idx-1) * exp(-dt(k_idx) / tau_discrete(i)) + ...
                                 gamma_est(i) * delta_theta * ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i)));
            end
        end
        % 시간 k_idx에서의 V_est 계산
        V_est(k_idx) = ocv_over_time(k_idx) + R0_est * ik(k_idx) + sum(V_RC(:, k_idx));
    end
    
    %% 4.5 DRT Gamma 그래프 출력
    figure(1);
    subplot(4, 4, s);
    plot(theta_discrete, gamma_est, 'LineWidth', 1.5);
    xlabel('\theta = ln(\tau)');
    ylabel('\gamma');
    title(['DRT for Trip ', num2str(s)]);
    grid on;
    
    %% 4.6 전압 비교 그래프 출력
    figure(2);
    subplot(4, 4, s);
    plot(t, V_sd, 'b', 'LineWidth', 1);
    hold on;
    plot(t, V_est, 'r--', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    title(['Voltage Comparison for Trip ', num2str(s)]);
    legend('Measured V_{sd}', 'Estimated V_{est}');
    grid on;
    hold off;
end


%% 5. Gamma(SOC, Theta) 3D 그래프 생성
% 각 트립의 SOC 중간값에 해당하는 Gamma 값을 3차원으로 배열
% soc_mid_all: (num_trips-1) x 1
% gamma_est_all: (num_trips-1) x n

% 정렬: SOC 중간값을 기준으로 오름차순 정렬
[soc_sorted, sort_idx] = sort(soc_mid_all);
gamma_sorted = gamma_est_all(sort_idx, :);

% 그리드 생성
[SOC_grid, Theta_grid] = meshgrid(soc_sorted, theta_discrete);

% Gamma 값을 전치하여 (n x num_trips-1) 행렬로 설정
Gamma_grid = gamma_sorted';

% 3D 서피스 플롯 생성
figure(3);
surf_handle = surf(SOC_grid, Theta_grid, Gamma_grid, ...
                  'FaceColor', 'blue', ...      % 단일 색상 지정
                  'EdgeColor', 'none');         % 테두리 색상 제거
xlabel('SOC');
ylabel('\theta = ln(\tau)');
zlabel('\gamma');
title('Gamma(SOC, \theta) 3D Surface Plot');
% colorbar;  % 색상 매핑이 없으므로 필요 없을 수 있음
view(135, 30);  % 시점 조정
grid on;

% 추가: 투명도 설정 (선택 사항)
alpha(0.8);

% 축의 방향을 맞추기 위해 X축과 Y축을 조정
axis tight;

%% save

save('gamma_data.mat', 'gamma_sorted', 'soc_sorted', 'theta_discrete', 'R0_est_all', 'soc_mid_all');
save('soc_ocv_data.mat', 'soc_values', 'ocv_values');


