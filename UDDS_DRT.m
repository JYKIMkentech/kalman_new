clear; clc; close all;

%% 1. UDDS 주행 데이터 로드
% UDDS 주행 데이터를 로드합니다.
load('C:\Users\deu04\OneDrive\바탕 화면\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
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
SOC_initial = 1;  % 초기 SOC를 100%로 가정합니다.

% 전류 데이터의 시간 간격(delta_t) 계산
delta_t = [0; diff(udds_time)];  % 각 측정 사이의 시간 간격

% 배터리 용량 (Ah를 Coulomb로 변환)
Q_battery = 2.9 * 3600;  % 배터리 용량 (2.9Ah)

% 시간에 따른 SOC 계산
udds_SOC = SOC_initial - cumsum(udds_current .* delta_t) / Q_battery;

% SOC 값이 0과 1 사이에 있도록 제한
udds_SOC(udds_SOC > 1) = 1;
udds_SOC(udds_SOC < 0) = 0;

%% 4. OCV 계산
% SOC에 기반하여 OCV 값을 보간합니다.
udds_OCV = interp1(soc_values, ocv_values, udds_SOC, 'linear', 'extrap');

%% 5. 파라미터 설정
n = 40;  % 이산 요소의 수
num_data = length(udds_time);  % 데이터 포인트의 수

% τ 값의 범위를 설정합니다 (0초에서 1e5초)
tau_min = 1e-5;  % 최소 τ 값 (0에 매우 가까운 값)
tau_max = 1e5;   % 최대 τ 값 (100,000초)

% τ 값을 로그 스케일로 생성합니다
tau_discrete = logspace(log10(tau_min), log10(tau_max), n);

% θ 값을 계산합니다
theta_discrete = log(tau_discrete);

% delta_theta 계산
delta_theta = theta_discrete(2) - theta_discrete(1);

%% 6. SOC 타겟 설정 및 γ와 R0 추정
% SOC 타겟을 5% 간격으로 설정합니다.
soc_targets = 1:-0.05:0;

% γ 값을 저장할 행렬을 NaN으로 초기화합니다.
gamma_matrix = NaN(length(soc_targets), n);

% R0 값을 저장할 벡터를 NaN으로 초기화합니다.
R0_vector = NaN(length(soc_targets), 1);

% 정규화 파라미터 λ 설정
lambda = 0.51795;  % 필요에 따라 조정

% 각 SOC 타겟에 대해 γ와 R0 추정
for soc_idx = 1:length(soc_targets)
    SOC_target = soc_targets(soc_idx);
    
    % SOC_target 주변의 데이터를 선택하기 위한 인덱스 찾기
    % SOC 변화가 너무 크지 않도록 범위를 조정합니다
    soc_tolerance = 0.02;  % 2%로 증가시킴
    indices = find(udds_SOC >= (SOC_target - soc_tolerance) & udds_SOC <= (SOC_target + soc_tolerance));
    
    % 충분한 데이터 포인트가 있는지 확인합니다.
    if length(indices) < 20  % 최소 데이터 포인트 수를 20개로 감소
        disp(['SOC ' num2str(SOC_target*100) '% 에서 데이터 포인트가 충분하지 않습니다.']);
        continue;
    end
    
    % 선택된 인덱스에서 데이터 추출
    idx_start = indices(1);
    idx_end = indices(end);
    
    current_window = udds_current(idx_start:idx_end);
    voltage_window = udds_voltage(idx_start:idx_end);
    time_window = udds_time(idx_start:idx_end);
    SOC_window = udds_SOC(idx_start:idx_end);
    OCV_window = udds_OCV(idx_start:idx_end);
    
    % delta_t 재계산
    delta_t_window = [0; diff(time_window)];
    
    % W 행렬 구성
    num_data_window = length(time_window);
    W_window = zeros(num_data_window, n);  % W 행렬 초기화
    
    for k = 1:num_data_window
        if k == 1
            for i = 1:n
                W_window(k, i) = current_window(k) * (1 - exp(-delta_t_window(k) / tau_discrete(i))) * delta_theta;
            end
        else
            for i = 1:n
                W_window(k, i) = W_window(k-1, i) * exp(-delta_t_window(k) / tau_discrete(i)) + ...
                                 current_window(k) * (1 - exp(-delta_t_window(k) / tau_discrete(i))) * delta_theta;
            end
        end
    end
    
    % y 벡터 계산
    y_window = voltage_window - OCV_window;
    
    % Φ 행렬 구성 (전류 및 W 행렬)
    Phi_window = [current_window, W_window];  % 전류를 첫 번째 열로 추가
    
    % 1차 차분 행렬 L 구성 (γ에만 적용)
    L = zeros(n-1, n);
    for i = 1:n-1
        L(i, i) = -1;
        L(i, i+1) = 1;
    end
    
    % 정규화 행렬 구성 (R0에 정규화를 적용하지 않음)
    L_extended = [zeros(n-1,1), L];  % L을 확장하여 R0 부분은 0으로 채웁니다.
    
    % 비용 함수 최소화를 위한 행렬 방정식 설정
    A = Phi_window' * Phi_window + lambda * (L_extended' * L_extended);
    b = Phi_window' * y_window;
    
    % θ 추정 (R0와 γ)
    theta_estimated = A \ b;
    
    % 추정된 R0와 γ 분리
    R0_estimated = theta_estimated(1);
    gamma_estimated = theta_estimated(2:end);
    
    % 결과 저장
    R0_vector(soc_idx) = R0_estimated;
    gamma_matrix(soc_idx, :) = gamma_estimated';
    
    disp(['SOC ' num2str(SOC_target*100) '% 에서 R0와 γ 추정 완료.']);
end

%% 7. 추정된 SOC 인덱스 찾기
valid_indices = ~isnan(R0_vector);

valid_soc_targets = soc_targets(valid_indices);
valid_gamma_matrix = gamma_matrix(valid_indices, :);
valid_R0_vector = R0_vector(valid_indices);

%% 8. γ(SOC, θ)의 3차원 그래프 그리기
% gamma_matrix의 행과 열에 해당하는 SOC와 θ 값을 생성합니다.
[Theta_grid, SOC_grid] = meshgrid(theta_discrete, valid_soc_targets);

% γ 값을 플로팅에 맞게 설정합니다.
gamma_matrix_plot = valid_gamma_matrix;

% 3D 그래프 그리기
figure;
surf(SOC_grid, Theta_grid, gamma_matrix_plot', 'EdgeColor', 'none');
xlabel('SOC');
ylabel('\theta = ln(\tau)');
zlabel('\gamma(\theta)');
title('γ(SOC, θ)의 3차원 그래프');
colorbar;
view(135, 30);  % 그래프 뷰 각도 조정
grid on;

%% 9. SOC에 따른 R0 값 플로팅
figure;
plot(valid_soc_targets, valid_R0_vector, 'o-', 'LineWidth', 1.5);
xlabel('SOC');
ylabel('R0 (Ω)');
title('SOC에 따른 R0 추정값');
grid on;

%% 10. 특정 SOC에서 측정된 전압과 모델 전압 비교
% 유효한 SOC 중 하나를 선택합니다 (예: 첫 번째 유효한 SOC)
soc_idx = 1;  % 또는 원하는 인덱스로 설정

if any(valid_gamma_matrix(soc_idx, :))
    % 해당 SOC에서의 데이터 인덱스 찾기
    SOC_target = valid_soc_targets(soc_idx);
    soc_tolerance = 0.02;  % 2%
    indices = find(udds_SOC >= (SOC_target - soc_tolerance) & udds_SOC <= (SOC_target + soc_tolerance));
    
    if ~isempty(indices)
        idx_start = indices(1);
        idx_end = indices(end);
        
        current_window = udds_current(idx_start:idx_end);
        voltage_window = udds_voltage(idx_start:idx_end);
        time_window = udds_time(idx_start:idx_end);
        SOC_window = udds_SOC(idx_start:idx_end);
        OCV_window = udds_OCV(idx_start:idx_end);
        
        % delta_t 재계산
        delta_t_window = [0; diff(time_window)];
        
        % W 행렬 재구성
        num_data_window = length(time_window);
        W_window = zeros(num_data_window, n);  % W 행렬 초기화
        
        for k = 1:num_data_window
            if k == 1
                for i = 1:n
                    W_window(k, i) = current_window(k) * (1 - exp(-delta_t_window(k) / tau_discrete(i))) * delta_theta;
                end
            else
                for i = 1:n
                    W_window(k, i) = W_window(k-1, i) * exp(-delta_t_window(k) / tau_discrete(i)) + ...
                                     current_window(k) * (1 - exp(-delta_t_window(k) / tau_discrete(i))) * delta_theta;
                end
            end
        end
        
        % 모델 전압 계산
        R0_estimated = valid_R0_vector(soc_idx);
        gamma_estimated = valid_gamma_matrix(soc_idx, :)';
        V_estimated = OCV_window + R0_estimated * current_window + W_window * gamma_estimated;
        
        % 측정된 전압과 모델 전압 비교
        figure;
        plot(time_window, voltage_window, 'k-', 'LineWidth', 1.5);
        hold on;
        plot(time_window, V_estimated, 'r--', 'LineWidth', 1.5);
        xlabel('시간 (s)');
        ylabel('전압 (V)');
        title(['측정 전압 vs. 모델 전압 (SOC ' num2str(SOC_target*100) '%)']);
        legend('측정 전압', '모델 전압');
        grid on;
    else
        disp(['SOC ' num2str(SOC_target*100) '% 근처에 데이터가 충분하지 않습니다.']);
    end
else
    disp(['SOC ' num2str(SOC_target*100) '%에서 γ 추정값이 없습니다.']);
end

