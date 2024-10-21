clc; clear; close all;

%% Parameters
n = 40;                      % RC 요소의 수
t = linspace(0, 100, 1000);  % 시간 벡터 (0초부터 100초까지, 1000포인트)
dt = t(2) - t(1);            % 시간 간격
num_scenarios = 10;          % 총 시나리오 수
validation_folds = {[1,2], [3,4], [5,6], [7,8]}; % 4개의 검증 세트
test_scenarios = [9,10];     % 테스트 세트
candidate_lambdas = logspace(-4, 4, 50); % \(\lambda\) 후보 (1e-4부터 1e4까지, 50개)

%% True DRT Parameters [theta, gamma]
mu_theta = log(10);  % 평균 ln(τ)
sigma_theta = 1;     % ln(τ)의 표준편차

% θ (ln(τ)) 값 설정
theta_min = mu_theta - 3*sigma_theta;
theta_max = mu_theta + 3*sigma_theta;
theta_discrete = linspace(theta_min, theta_max, n);  % ln(τ)를 균등 간격으로 분할

% τ 값은 θ의 지수 함수
tau_discrete = exp(theta_discrete);  % τ = exp(θ)

% Δθ 계산
delta_theta = theta_discrete(2) - theta_discrete(1);  % Δθ = θ_{n+1} - θ_n

%% Regularization Matrix L (First-order difference, no scaling)
% 차분 행렬 D 생성
D = zeros(n-1, n);
for i = 1:n-1
    D(i, i) = -1;
    D(i, i+1) = 1;
end
L = D;  % 스케일링 적용하지 않음

%% True DRT [theta, gamma]
g_discrete_true = normpdf(theta_discrete, mu_theta, sigma_theta);
gamma_discrete_true = g_discrete_true;  % γ(theta) = g(theta)
gamma_discrete_true = gamma_discrete_true / max(gamma_discrete_true);  % 최대값으로 정규화

%% Define Amplitudes and Periods for Current Synthesis
A = [1, 1, 1;          % 시나리오 1
     1.7, 0.6, 0.7;    % 시나리오 2
     0.2, 0.5, 2.3;    % 시나리오 3
     1.3, 1.1, 0.6;    % 시나리오 4
     1.7, 1.8, 0.5;    % 시나리오 5
     1.27, 1.33, 0.4;  % 시나리오 6
     1.2, 1.6, 0.2;    % 시나리오 7
     0.9, 0.7, 2.4;    % 시나리오 8
     1.1, 1.1, 0.8;    % 시나리오 9
     0.1, 0.1, 2.8];   % 시나리오 10

T = [1, 5, 20;         % 시나리오 1
     2, 4, 20;         % 시나리오 2
     1, 20, 25;        % 시나리오 3
     1.5, 5.3, 19.8;   % 시나리오 4
     2.5, 4.2, 20.5;   % 시나리오 5
     1.5, 20.9, 24.2;  % 시나리오 6
     1.3, 6, 19.3;     % 시나리오 7
     2.2, 4.8, 20.2;   % 시나리오 8
     2, 20.8, 26.1;    % 시나리오 9
     1.1, 4.3, 20.1];  % 시나리오 10

%% Generate Synthetic Current Data (Multi-Sine Approach)
ik_scenarios = zeros(num_scenarios, length(t)); % 시나리오별 전류 초기화

for s = 1:num_scenarios
    % 각 시나리오에 대해 세 개의 사인파 합성
    ik_scenarios(s, :) = A(s,1)*sin(2*pi*t / T(s,1)) + ...
                         A(s,2)*sin(2*pi*t / T(s,2)) + ...
                         A(s,3)*sin(2*pi*t / T(s,3));
end

%% Generate Voltage Data for Each Scenario
V_est_all = zeros(num_scenarios, length(t)); % 추정 전압 저장
V_sd_all = cell(num_scenarios, 1);          % 노이즈가 추가된 전압 저장

rng(0); % 동일한 노이즈 패턴 적용
noise_level = 0.01; % 노이즈 수준

for s = 1:num_scenarios
    fprintf('시나리오 %d/%d 데이터 생성 중...\n', s, num_scenarios);
    
    ik = ik_scenarios(s, :); % 현재 시나리오의 전류
    
    %% 전압 초기화
    V_est = zeros(1, length(t)); % 추정 전압
    R0 = 0.1; % 오믹 저항
    OCV = 0;  % 개방 회로 전압
    V_RC = zeros(n, length(t)); % 각 RC 요소의 전압
    
    %% 초기 전압 계산 (k=1)
    for i = 1:n
        V_RC(i, 1) = gamma_discrete_true(i) * delta_theta * ik(1) * (1 - exp(-dt / tau_discrete(i)));
    end
    V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));
    
    %% 시간 t > 1에 대한 전압 계산
    for k_idx = 2:length(t)
        for i = 1:n
            V_RC(i, k_idx) = V_RC(i, k_idx-1) * exp(-dt / tau_discrete(i)) + ...
                              gamma_discrete_true(i) * delta_theta * ik(k_idx) * (1 - exp(-dt / tau_discrete(i)));
        end
        V_est(k_idx) = OCV + R0 * ik(k_idx) + sum(V_RC(:, k_idx));
    end
    
    % 추정 전압 저장 및 노이즈 추가
    V_est_all(s, :) = V_est;
    V_sd = V_est + noise_level * randn(size(V_est));
    V_sd_all{s} = V_sd;
    
    %% 시나리오별 W 행렬 구성
    W = zeros(length(t), n); % W 행렬 초기화
    for k_idx = 1:length(t)
        for i = 1:n
            if k_idx == 1
                W(k_idx, i) = ik(k_idx) * (1 - exp(-dt / tau_discrete(i))) * delta_theta;
            else
                W(k_idx, i) = W(k_idx-1, i) * exp(-dt / tau_discrete(i)) + ...
                              ik(k_idx) * (1 - exp(-dt / tau_discrete(i))) * delta_theta;
            end
        end
    end
    W_all{s} = W; % W 행렬 저장
end

%% Cross-Validation to Optimize Lambda and Estimate Uncertainty
CVE = zeros(length(candidate_lambdas),1); % 각 \(\lambda\)에 대한 교차 검증 오류 초기화
gamma_fold_collection = cell(length(candidate_lambdas), length(validation_folds)); % 각 폴드에서의 gamma 저장

for l = 1:length(candidate_lambdas)
    current_lambda = candidate_lambdas(l); % 현재 \(\lambda\)
    total_error = 0; % 총 교차 검증 오류 초기화
    
    for fold = 1:length(validation_folds)
        val_set = validation_folds{fold}; % 현재 검증 세트
        train_set = setdiff(1:num_scenarios, [val_set, test_scenarios]); % 훈련 세트 (테스트 세트 제외)
        
        % 훈련 세트에서 W_train과 y_train 구성
        W_train = [];
        y_train = [];
        
        for s = train_set
            W_train = [W_train; W_all{s}];
            y_train = [y_train; V_sd_all{s}' - 0.1 * ik_scenarios(s, :)' - 0]; % OCV=0
        end
        
        %% 정규화된 최소 제곱 해 계산
        gamma_analytical = (W_train' * W_train + current_lambda * (L' * L)) \ (W_train' * y_train);
        gamma_analytical(gamma_analytical < 0) = 0; % 음수 값 제거
        
        % 각 폴드에서의 gamma 저장
        gamma_fold_collection{l, fold} = gamma_analytical;
        
        %% 검증 세트에서 오류 계산
        for s_val = val_set
            W_val = W_all{s_val};
            V_val = V_sd_all{s_val};
            ik_val = ik_scenarios(s_val, :)';
            y_val = V_val' - 0.1 * ik_val - 0; % OCV=0
            
            % 예측 전압 계산
            V_pred = W_val * gamma_analytical + 0.1 * ik_val + 0;
            
            % 제곱 오차 합산
            error = sum((V_val' - V_pred).^2);
            total_error = total_error + error;
        end
    end
    
    % 현재 \(\lambda\)에 대한 CVE 저장
    CVE(l) = total_error;
end

%% Find Optimal Lambda
[~, optimal_idx] = min(CVE);
optimal_lambda = candidate_lambdas(optimal_idx);
fprintf('최적 \lambda: %e\n', optimal_lambda);

%% Compute Uncertainty from Cross-Validation
% 최적 \(\lambda\)에 대한 각 폴드의 gamma 수집
gamma_folds = gamma_fold_collection(optimal_idx, :)';
gamma_matrix = zeros(length(gamma_folds), n); % (폴드 수 x n) 크기의 행렬

for k = 1:length(gamma_folds)
    gamma_vector = gamma_folds{k};
    if size(gamma_vector, 2) == 1
        gamma_vector = gamma_vector'; % 행 벡터로 변환
    end
    gamma_matrix(k, :) = gamma_vector;
end

% gamma의 평균과 표준편차 계산
gamma_mean = mean(gamma_matrix, 1); % (1 x n)
gamma_std = std(gamma_matrix, 0, 1); % (1 x n)

%% Plot CVE vs Lambda with Logarithmic Scale and Appropriate Y-axis Limits
figure;
semilogx(candidate_lambdas, CVE, 'b-', 'LineWidth', 1.5); % CVE vs Lambda plot
hold on;

% 최적 \(\lambda\) 포인트 표시
semilogx(optimal_lambda, CVE(optimal_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% 최적 \(\lambda\) 텍스트 추가
optimal_lambda_str = ['최적 \lambda = ', num2str(optimal_lambda, '%.2e')];

% 레이블 및 제목 (이스케이프 문자 수정)
xlabel('\\lambda (정규화 파라미터)');
ylabel('교차 검증 오류 (CVE)');
title('로그 스케일 \lambda 에 따른 CVE 그래프');

% 그리드 및 범례
grid on;
set(gca, 'YScale', 'log');  % Y축 로그 스케일 설정
% Y축 한계 조정 (데이터 범위에 맞게 설정)
ylim([min(CVE(CVE > 0)) * 0.9, max(CVE) * 1.1]);
legend({'CVE', optimal_lambda_str}, 'Location', 'best');
hold off;

%% Plot DRT with Uncertainty (Error Bars)
figure;
errorbar(theta_discrete, gamma_mean, gamma_std, 'r--o', 'LineWidth', 1.5, 'DisplayName', 'Estimated \gamma(\theta) with Uncertainty');
hold on;
plot(theta_discrete, gamma_discrete_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True \gamma(\theta)');
xlabel('ln(\tau) = \theta');
ylabel('\gamma(\theta)');
title('True DRT와 추정된 DRT (불확실성 포함) 비교');
legend('Location', 'Best');
grid on;
hold off;

% 그래프 레이아웃 정규화 설정
set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);  % 화면 비율로 설정

%% Retrain on All Training Data with Optimal Lambda and Evaluate on Test Set
% 훈련 세트: 시나리오 [1-8]
train_set_final = setdiff(1:num_scenarios, test_scenarios); % [1-8]
sum_WtW_final = zeros(n, n);   % W^T * W 합 초기화
sum_WtV_final = zeros(n, 1);   % W^T * y 합 초기화

for s = train_set_final
    W_s = W_all{s};
    V_s = V_sd_all{s};
    ik_s = ik_scenarios(s, :)';
    y_s = V_s' - 0.1 * ik_s - 0;  % OCV=0
    sum_WtW_final = sum_WtW_final + W_s' * W_s;
    sum_WtV_final = sum_WtV_final + W_s' * y_s;
end

% 정규화 항
regularization_term = optimal_lambda * (L' * L);

% 최적 \(\gamma\) 계산
gamma_final = (sum_WtW_final + regularization_term) \ sum_WtV_final;
gamma_final(gamma_final < 0) = 0; % 음수 값 제거

%% Test on Test Scenarios [9,10]
test_set = test_scenarios;
total_test_error = 0;

for s_test = test_set
    W_test = W_all{s_test};
    V_test = V_sd_all{s_test};
    ik_test = ik_scenarios(s_test, :)';
    y_test = V_test' - 0.1 * ik_test - 0; % OCV=0
    
    % 예측 전압 계산
    V_pred_test = W_test * gamma_final + 0.1 * ik_test + 0;
    
    % 제곱 오차 합산
    error_test = sum((V_test' - V_pred_test).^2);
    total_test_error = total_test_error + error_test;
end

fprintf('테스트 세트에서의 총 제곱 전압 오류: %f\n', total_test_error);

%% Plot DRT Comparison: True vs Optimized Gamma (Without Uncertainty)
figure;
plot(theta_discrete, gamma_discrete_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True \gamma(\theta)');
hold on;
plot(theta_discrete, gamma_final, 'b--', 'LineWidth', 2, 'DisplayName', 'Optimized \gamma(\theta)');
xlabel('ln(\tau) = \theta');
ylabel('\gamma(\theta)');
title('True DRT와 최적화된 DRT 비교 (불확실성 미포함)');
legend('Location', 'Best');
grid on;
hold off;

% 그래프 레이아웃 정규화 설정
set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);  % 화면 비율로 설정


