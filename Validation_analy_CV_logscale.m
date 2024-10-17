clc; clear; close all;

%% Parameters
n = 21;                      % Number of RC elements
t = linspace(0, 100, 1000);  % Time vector (0 to 100 seconds with 1000 points)
dt = t(2) - t(1);            % Time step
num_scenarios = 10;          % Total number of scenarios
test_scenarios = [9, 10];    % Test scenarios
validation_folds = {[1,2], [3,4], [5,6], [7,8]}; % 4 validation sets
candidate_lambdas = logspace(-4, 4, 50); % Lambda candidates (1e-4 to 1e4)

%% True DRT Parameters [theta, gamma]
mu_theta = log(10);  % 평균 ln(τ)
sigma_theta = 0.2;   % ln(τ)의 표준편차

% θ (ln(τ)) 값 설정
theta_min = mu_theta - 3*sigma_theta;
theta_max = mu_theta + 3*sigma_theta;
theta_discrete = linspace(theta_min, theta_max, n);  % ln(τ) 값을 일정 간격으로 분할

% τ 값은 θ의 지수 함수
tau_discrete = exp(theta_discrete);  % τ = exp(θ)

% Δτ_log 계산
delta_tau_log = theta_discrete(2) - theta_discrete(1);  % Δτ_log = ln(τ_{n+1}) - ln(τ_n)

%% Regularization Matrix L (First-order derivative)
% Difference matrix D 생성
D = zeros(n-1, n);
for i = 1:n-1
    D(i, i) = -1;
    D(i, i+1) = 1;
end

% 정규화 행렬 L 생성 (스케일링 적용)
L = (1 / delta_tau_log) * D;

%% True DRT [θ, gamma]
g_discrete_true = normpdf(theta_discrete, mu_theta, sigma_theta);
gamma_discrete_true = g_discrete_true;  % γ(θ) = g(τ)
gamma_discrete_true = gamma_discrete_true / max(gamma_discrete_true);  % 최대값으로 정규화

% R_n 계산
R_n_true = gamma_discrete_true * delta_tau_log;  % R_n = γ_n * Δτ_log

%% Initialize Storage Variables
W_all = cell(num_scenarios, 1);      % Store W matrices for all scenarios
V_sd_all = cell(num_scenarios, 1);   % Store noisy V data for all scenarios

%% Define Amplitudes and Periods for Current Synthesis
A = [1, 1, 1;          % Scenario 1
     1.7, 0.6, 0.7;    % Scenario 2
     0.2, 0.5, 2.3;    % Scenario 3
     1.3, 1.1, 0.6;    % Scenario 4
     1.7, 1.8, 0.5;    % Scenario 5
     1.27, 1.33, 0.4;  % Scenario 6
     1.2, 1.6, 0.2;    % Scenario 7
     0.9, 0.7, 2.4;    % Scenario 8
     1.1, 1.1, 0.8;    % Scenario 9
     0.1, 0.1, 2.8];   % Scenario 10

T = [1, 5, 20;         % Scenario 1
     2, 4, 20;         % Scenario 2
     1, 20, 25;        % Scenario 3
     1.5, 5.3, 19.8;   % Scenario 4
     2.5, 4.2, 20.5;   % Scenario 5
     1.5, 20.9, 24.2;  % Scenario 6
     1.3, 6, 19.3;     % Scenario 7
     2.2, 4.8, 20.2;   % Scenario 8
     2, 20.8, 26.1;    % Scenario 9
     1.1, 4.3, 20.1];  % Scenario 10

%% Generate Synthetic Current Data (Multi-Sine Approach)
ik_scenarios = zeros(num_scenarios, length(t)); % Initialize current scenarios

for s = 1:num_scenarios
    % Sum of three sine waves for each scenario
    ik_scenarios(s, :) = A(s,1)*sin(2*pi*t / T(s,1)) + ...
                         A(s,2)*sin(2*pi*t / T(s,2)) + ...
                         A(s,3)*sin(2*pi*t / T(s,3));
end

%% Generate Voltage Data for Each Scenario
V_est_all = zeros(num_scenarios, length(t)); % Store estimated V
rng(0); % For reproducibility
noise_level = 0.01; % Noise level

for s = 1:num_scenarios
    fprintf('Generating data for Scenario %d/%d...\n', s, num_scenarios);
    
    ik = ik_scenarios(s, :); % Current for scenario s
    
    %% Initialize Voltage
    V_est = zeros(1, length(t)); % Estimated voltage
    R0 = 0.1; % Ohmic resistance
    OCV = 0;  % Open Circuit Voltage
    V_RC = zeros(n, length(t)); % RC voltages for each element
    
    %% Initial Voltage Calculation (k=1)
    for i = 1:n
        V_RC(i, 1) = ik(1) * R_n_true(i) * (1 - exp(-dt / tau_discrete(i)));
    end
    V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));
    
    %% Voltage Calculation for t > 1
    for k_idx = 2:length(t)
        for i = 1:n
            V_RC(i, k_idx) = exp(-dt / tau_discrete(i)) * V_RC(i, k_idx-1) + ...
                              R_n_true(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k_idx);
        end
        V_est(k_idx) = OCV + R0 * ik(k_idx) + sum(V_RC(:, k_idx));
    end
    
    % Store V_est and add noise
    V_est_all(s, :) = V_est;
    V_sd = V_est + noise_level * randn(size(V_est));
    V_sd_all{s} = V_sd;
    
    %% Construct W Matrix for Scenario s
    W = zeros(length(t), n); % Initialize W matrix
    for k_idx = 1:length(t)
        for i = 1:n
            if k_idx == 1
                W(k_idx, i) = (1 - exp(-dt / tau_discrete(i))) * ik(k_idx);
            else
                W(k_idx, i) = exp(-dt / tau_discrete(i)) * W(k_idx-1, i) + ...
                              (1 - exp(-dt / tau_discrete(i))) * ik(k_idx);
            end
        end
    end
    % W 행렬에 Δτ_log를 곱해줌
    W = W * delta_tau_log;
    W_all{s} = W; % Store W matrix
end

%% Cross-Validation to Optimize Lambda
CVE = zeros(length(candidate_lambdas),1); % Initialize CVE for each lambda

for l = 1:length(candidate_lambdas)
    current_lambda = candidate_lambdas(l); % Current lambda
    total_error = 0; % Initialize total CVE
    
    for fold = 1:length(validation_folds)
        val_set = validation_folds{fold}; % Current validation set
        
        % Define training set: exclude test and current validation set
        train_set = setdiff(1:num_scenarios, [test_scenarios, val_set]);
        
        % Stack W_train and V_train from training set
        sum_WtW = zeros(n, n);   % Initialize sum(W_s^T * W_s)
        sum_WtV = zeros(n, 1);   % Initialize sum(W_s^T * V_s)
        
        for s = train_set
            W_s = W_all{s};
            V_s = V_sd_all{s};
            ik_s = ik_scenarios(s, :)';
            y_s = V_s' - R0 * ik_s - OCV;  % Adjusted voltage
            sum_WtW = sum_WtW + W_s' * W_s;
            sum_WtV = sum_WtV + W_s' * y_s;
        end
        
        %% Regularization term with scaling
        regularization_term = current_lambda * (L' * L);
        
        %% Analytical Solution for R_train
        R_train = (sum_WtW + regularization_term) \ sum_WtV;
        R_train(R_train < 0) = 0; % Enforce non-negativity
        
        %% Validation Error Calculation
        for s_val = val_set
            W_val = W_all{s_val};
            V_val = V_sd_all{s_val};
            ik_val = ik_scenarios(s_val, :)';
            y_val = V_val' - R0 * ik_val - OCV;  % Adjusted voltage
            
            % Predict V_pred for validation set
            V_pred = W_val * R_train + R0 * ik_val + OCV;
            
            % Compute squared voltage error
            error = sum((V_val' - V_pred).^2);
            total_error = total_error + error;
        end
    end
    
    % Store CVE for current lambda
    CVE(l) = total_error;
end

%% Find Optimal Lambda
[~, optimal_idx] = min(CVE);
optimal_lambda = candidate_lambdas(optimal_idx);
fprintf('Optimal lambda: %e\n', optimal_lambda);

%% Plot CVE vs Lambda with Aggressive Y-axis Zoom and Logarithmic Scaling
figure;
semilogx(candidate_lambdas, CVE, 'b-', 'LineWidth', 1.5); % CVE vs Lambda plot
hold on;

% Plot optimal lambda point
semilogx(optimal_lambda, CVE(optimal_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Optimal lambda point

% Add text to indicate the optimal lambda inside the legend
optimal_lambda_str = ['Optimal \lambda = ', num2str(optimal_lambda, '%.2e')]; % Text string for optimal lambda

% Labels and title
xlabel('\lambda');
ylabel('Crossvalidation Error (CVE)');
title('CVE vs \lambda');

% Grid and legend
grid on;

% Set Y-axis to logarithmic scale
set(gca, 'YScale', 'log');  % Set Y-axis to log scale

% Adjust Y-axis limits
ylim([min(CVE(CVE > 0)) * 0.9, max(CVE) * 1.1]);

% Add legend and include optimal lambda in the legend
legend({'CVE', optimal_lambda_str}, 'Location', 'best');

hold off;

%% Retrain on All Training Data with Optimal Lambda
% Training set: exclude test scenarios (9 and 10)
train_set_final = setdiff(1:num_scenarios, test_scenarios);
sum_WtW_final = zeros(n, n);   % Initialize sum(W_s^T * W_s)
sum_WtV_final = zeros(n, 1);   % Initialize sum(W_s^T * V_s)

for s = train_set_final
    W_s = W_all{s};
    V_s = V_sd_all{s};
    ik_s = ik_scenarios(s, :)';
    y_s = V_s' - R0 * ik_s - OCV;  % Adjusted voltage
    sum_WtW_final = sum_WtW_final + W_s' * W_s;
    sum_WtV_final = sum_WtV_final + W_s' * y_s;
end

% Regularization term with scaling
regularization_term = optimal_lambda * (L' * L);

% Compute R_final using optimal lambda
R_final = (sum_WtW_final + regularization_term) \ sum_WtV_final;
R_final(R_final < 0) = 0; % Enforce non-negativity

% γ_n 추정
gamma_final = R_final / delta_tau_log;

%% Test on Test Scenarios (9 and 10)
test_set = test_scenarios;
total_test_error = 0;

for s_test = test_set
    W_test = W_all{s_test};
    V_test = V_sd_all{s_test};
    ik_test = ik_scenarios(s_test, :)';
    y_test = V_test' - R0 * ik_test - OCV;  % Adjusted voltage
    
    % Predict V_pred_test
    V_pred_test = W_test * R_final + R0 * ik_test + OCV;
    
    % Compute squared voltage error
    error_test = sum((V_test' - V_pred_test).^2);
    total_test_error = total_test_error + error_test;
end

fprintf('Total squared voltage error on test data: %f\n', total_test_error);

%% Plot DRT Comparison: True vs Optimized Gamma
figure;
plot(theta_discrete, gamma_discrete_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True \gamma(\theta)');
hold on;
plot(theta_discrete, gamma_final, 'r--', 'LineWidth', 2, 'DisplayName', 'Optimized \gamma(\theta)');
xlabel('ln(\tau) = \theta');
ylabel('\gamma(\theta)');
title('Comparison of True DRT and Optimized DRT');
legend('Location', 'Best');
grid on;
hold off;

% 그래프 레이아웃 정규화 설정
set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);  % 화면 비율로 설정

