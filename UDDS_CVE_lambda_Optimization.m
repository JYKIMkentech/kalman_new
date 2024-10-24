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

% 정규화 파라미터 람다 값 범위 설정
lambda_values = logspace(-4, 3 , 10);  % 필요에 따라 조정 가능

% Gamma에 대한 1차 차분 행렬 L_gamma 생성
L_gamma = zeros(n-1, n);
for i = 1:n-1
    L_gamma(i, i) = -1;
    L_gamma(i, i+1) = 1;
end

% R0에 대한 정규화를 피하기 위해 L_aug 생성
L_aug = [L_gamma, zeros(n-1, 1)];

%% 4. 각 사이클의 데이터 준비
num_cycles = length(udds_data)-3;

% 사이클별 데이터 저장을 위한 셀 배열 초기화
t_all = cell(num_cycles, 1);
ik_all = cell(num_cycles, 1);
V_sd_all = cell(num_cycles, 1);
SOC_all = cell(num_cycles, 1);

for s = 1:num_cycles
    % 현재 사이클의 데이터 추출
    t = udds_data(s).t;
    ik = udds_data(s).I;
    V_sd = udds_data(s).V;
    SOC = udds_data(s).SOC;
    
    % 데이터 저장
    t_all{s} = t;
    ik_all{s} = ik;
    V_sd_all{s} = V_sd;
    SOC_all{s} = SOC;
end

%% 5. 교차 검증을 통한 람다 최적화

% 전체 사이클에서 2개를 검증 세트로 선택하는 조합 생성
validation_indices = nchoosek(1:num_cycles, 2);  % 16C2 = 120개의 조합
num_folds = size(validation_indices, 1);  % 120개의 폴드

cve_lambda = zeros(length(lambda_values), 1);  % 각 람다에 대한 CVE 저장

for l_idx = 1:length(lambda_values)
    lambda = lambda_values(l_idx);
    cve_total = 0;
    
    fprintf('Processing lambda %e (%d/%d)...\n', lambda, l_idx, length(lambda_values));
    
    for fold = 1:num_folds
        % 검증 세트와 학습 세트 분리
        val_scenarios = validation_indices(fold, :);
        train_scenarios = setdiff(1:num_cycles, val_scenarios);
        
        % 학습 데이터로 gamma 및 R0 추정
        [gamma_estimated, R0_estimated] = estimate_gamma(lambda, train_scenarios, ik_all, V_sd_all, SOC_all, soc_values, ocv_values, tau_discrete, delta_theta, L_aug);
        
        % 검증 데이터로 전압 예측 및 에러 계산
        error_fold = calculate_error(gamma_estimated, R0_estimated, val_scenarios, ik_all, V_sd_all, SOC_all, soc_values, ocv_values, tau_discrete, delta_theta);
        
        % 폴드의 에러 합산
        cve_total = cve_total + error_fold;
    end
    
    % 평균 CVE 계산
    cve_lambda(l_idx) = cve_total / num_folds ;
    fprintf('Lambda %e, CVE: %f\n', lambda, cve_lambda(l_idx));
end

%% CVE vs 람다 그래프 그리기
figure;
semilogx(lambda_values, cve_lambda, 'b-', 'LineWidth', 1.5); % CVE vs Lambda plot
hold on;

% 최적 \(\lambda\) 포인트 찾기
[~, min_idx] = min(cve_lambda);
optimal_lambda = lambda_values(min_idx);

% 최적 \(\lambda\) 포인트 표시
semilogx(optimal_lambda, cve_lambda(min_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% 최적 \(\lambda\) 텍스트 추가
optimal_lambda_str = ['Optimal \lambda = ', num2str(optimal_lambda, '%.2e')];
%ylim([19.7263, 19.7275]);

% 레이블 및 제목
xlabel('\lambda (Regularization Parameter)');
ylabel('Cross-Validation Error (CVE)');
title('CVE vs. \lambda (Log Scale)');

% 그리드 및 범례
grid on;
legend({'CVE', optimal_lambda_str}, 'Location', 'best');

hold off;

%% 최적의 람다로 전체 데이터로 gamma 및 R0 추정
[gamma_optimal, R0_optimal] = estimate_gamma(optimal_lambda, 1:num_cycles, ik_all, V_sd_all, SOC_all, soc_values, ocv_values, tau_discrete, delta_theta, L_aug);

%% 함수 정의

% Gamma 및 R0 추정 함수
function [gamma_estimated, R0_estimated] = estimate_gamma(lambda, train_scenarios, ik_all, V_sd_all, SOC_all, soc_values, ocv_values, tau_discrete, delta_theta, L_aug)
    n = length(tau_discrete);
    num_train = length(train_scenarios);
    
    gamma_estimated = zeros(n, num_train);
    R0_estimated = zeros(num_train, 1);
    
    for idx = 1:num_train
        s = train_scenarios(idx);
        ik = ik_all{s};
        V_sd = V_sd_all{s};
        SOC = SOC_all{s};
        
        % 시간 벡터
        t = (0:length(ik)-1)';  % 가정: dt가 일정하거나 필요에 따라 수정
        
        % 시간 간격 dt 계산
        delta_t = [0; diff(t)];
        dt = delta_t;
        if dt(1) == 0  % 첫 번째 dt 값이 0이면
            dt(1) = dt(2);  % 두 번째 dt 값으로 대체
        end
        
        % OCV 계산
        ocv_over_time = interp1(soc_values, ocv_values, SOC, 'linear', 'extrap');
        
        % W 행렬 생성
        W = compute_W(ik, tau_discrete, delta_theta, dt);
        
        % R0 추정을 위한 W 행렬 확장
        W_aug = [W, ik(:)];  % ik(:)는 ik를 열 벡터로 변환
        
        % y 벡터 생성
        y = V_sd - ocv_over_time;
        y = y(:);  % y를 열 벡터로 변환
        
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
            warning('Optimization did not converge for cycle %d.', s);
        end
        
        % gamma와 R0 추정값 추출
        gamma_est = Theta_est(1:n);
        R0_est = Theta_est(n+1);
        
        % 추정값 저장
        gamma_estimated(:, idx) = gamma_est;
        R0_estimated(idx) = R0_est;
    end
end

% 에러 계산 함수
function error_total = calculate_error(gamma_estimated, R0_estimated, val_scenarios, ik_all, V_sd_all, SOC_all, soc_values, ocv_values, tau_discrete, delta_theta)
    error_total = 0;
    num_val = length(val_scenarios);
    
    for idx = 1:num_val
        s = val_scenarios(idx);
        ik = ik_all{s};
        V_sd = V_sd_all{s};
        SOC = SOC_all{s};
        
        % OCV 계산
        %ocv_over_time = interp1(soc_values, ocv_values, SOC, 'linear', 'extrap');
        
        % 추정된 전압 계산
        gamma_est = mean(gamma_estimated, 2);  % 학습 세트의 평균 gamma 사용
        R0_est = mean(R0_estimated);          % 학습 세트의 평균 R0 사용
        V_est = predict_voltage(gamma_est, R0_est, ik, SOC, soc_values, ocv_values, tau_discrete, delta_theta);
        
        % V_sd와 V_est를 열 벡터로 변환 (필요한 경우)
        V_sd = V_sd(:);
        V_est = V_est(:);
        
        % 전압 차이의 제곱 합산
        error_total = error_total + sum((V_est - V_sd).^2);
    end
end

% W 행렬 계산 함수
function W = compute_W(ik, tau_discrete, delta_theta, dt)
    n = length(tau_discrete);
    len_t = length(ik);
    W = zeros(len_t, n);  % Initialize W matrix
    
    for k_idx = 1:len_t
        if k_idx == 1
            for i = 1:n
                W(k_idx, i) = ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i))) * delta_theta;
            end
        else
            for i = 1:n
                W(k_idx, i) = W(k_idx-1, i) * exp(-dt(k_idx) / tau_discrete(i)) + ...
                              ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i))) * delta_theta;
            end
        end
    end
end

% 전압 예측 함수
function V_est = predict_voltage(gamma_est, R0_est, ik, SOC, soc_values, ocv_values, tau_discrete, delta_theta)
    len_t = length(ik);
    n = length(gamma_est);
    V_RC = zeros(n, len_t);
    V_est = zeros(len_t, 1);
    
    % 시간 벡터
    t = (0:len_t-1)';  % 가정: dt가 일정하거나 필요에 따라 수정
    
    % 시간 간격 dt 계산
    delta_t = [0; diff(t)];
    dt = delta_t;
    if dt(1) == 0  % 첫 번째 dt 값이 0이면
        dt(1) = dt(2);  % 두 번째 dt 값으로 대체
    end
    
    % OCV 계산
    ocv_over_time = interp1(soc_values, ocv_values, SOC, 'linear', 'extrap');
    
    for k_idx = 1:len_t
        if k_idx == 1
            for i = 1:n
                V_RC(i, k_idx) = gamma_est(i) * ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i))) * delta_theta;
            end
        else
            for i = 1:n
                V_RC(i, k_idx) = V_RC(i, k_idx-1) * exp(-dt(k_idx) / tau_discrete(i)) + ...
                                 gamma_est(i) * ik(k_idx) * (1 - exp(-dt(k_idx) / tau_discrete(i))) * delta_theta;
            end
        end
        V_est(k_idx) = ocv_over_time(k_idx) + R0_est * ik(k_idx) + sum(V_RC(:, k_idx));
    end
end 