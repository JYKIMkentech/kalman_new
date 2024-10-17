clc; clear; close all;

%% Parameters 
n = 21;  % Number of RC elements
t = 0:0.01:100;  % Time vector (0 to 100 seconds with 10001 points)
dt = t(2) - t(1); % Time step
num_scenarios = 10;  % Number of current scenarios
lambda = 0.05;  % Regularization parameter

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

%% DRT Parameters

% ln(τ) = θ가 정규분포를 따름
mu_theta = log(10);  % ln(τ)의 평균
sigma_theta = 0.2;   % ln(τ)의 표준편차

% θ (ln(τ)) 값 설정
theta_min = mu_theta - 3*sigma_theta;
theta_max = mu_theta + 3*sigma_theta;
theta_discrete = linspace(theta_min, theta_max, n);  % ln(τ) 값을 일정 간격으로 분할

% τ 값은 θ의 지수 함수
tau_discrete = exp(theta_discrete);  % τ = exp(θ)

% Δτ_log 계산
delta_tau_log = theta_discrete(2) - theta_discrete(1);  % Δτ_log = ln(τ_{n+1}) - ln(τ_n)

% 참 DRT [θ, γ]
g_discrete_true = normpdf(theta_discrete, mu_theta, sigma_theta);
gamma_discrete_true = g_discrete_true;  % γ(θ) = g(τ)
gamma_discrete_true = gamma_discrete_true / max(gamma_discrete_true);  % 최대값으로 정규화

% R_n 계산
R_n_true = gamma_discrete_true * delta_tau_log;  % R_n = γ_n * Δτ_log

%% Regularization Matrix L (First-order derivative)
% Difference matrix D 생성
D = zeros(n-1, n);
for i = 1:n-1
    D(i, i) = -1;
    D(i, i+1) = 1;
end

% 정규화 행렬 L 생성 (스케일링 적용)
L = (1 / delta_tau_log) * D;

%% Initialize Storage Variables
W_all = cell(num_scenarios, 1);      % Store W matrices for all scenarios
V_sd_all = zeros(num_scenarios, length(t));   % Store noisy V data for all scenarios

%% Plot Synthesized Currents and Voltages
figure;
sgtitle('Synthesized Current and Voltage for Each Scenario');

% Initialize subplots with current plots
for s = 1:num_scenarios
    subplot(5, 2, s);
    yyaxis left
    plot(t, ik_scenarios(s, :), 'b-', 'LineWidth', 1.5);
    ylabel('Current (A)');
    xlabel('Time (s)');
    grid on;
end

%% Voltage Synthesis
for s = 1:num_scenarios
    fprintf('Processing Scenario %d/%d...\n', s, num_scenarios);
    
    % Current for the scenario
    ik = ik_scenarios(s, :);  % Current scenario input
    
    %% Initialize Voltage
    V_est = zeros(1, length(t));  % Model voltage calculated via n-RC model
    R0 = 0.1;  % Ohmic resistance
    OCV = 0;   % Open Circuit Voltage
    V_RC = zeros(n, length(t));  % RC voltages for each element
    
    %% Initial Voltage Calculation (first time step)
    for i = 1:n
        V_RC(i, 1) = R_n_true(i) * (1 - exp(-dt / tau_discrete(i))) * ik(1);
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
    
    % Store V_est for the current scenario
    V_est_all(s, :) = V_est;  % Save the calculated V_est for this scenario
    
    %% Add Noise to the Voltage
    rng(s);  % Ensure reproducibility of noise for each scenario
    noise_level = 0.01;
    V_sd = V_est + noise_level * randn(size(V_est));  % V_sd = synthetic measured voltage
    
    % Store V_sd for the current scenario
    V_sd_all(s, :) = V_sd;  % Save the noisy V_sd for this scenario
    
    %% Plot Voltage on Existing Subplots
    subplot(5, 2, s);
    yyaxis right
    plot(t, V_sd, 'r-', 'LineWidth', 1.5);
    ylabel('Voltage (V)');
    ylim([min(V_sd)-0.1, max(V_sd)+0.1]);
    
    % Update title with correct amplitudes and periods
    title(['Scenario ', num2str(s), ...
           ': A1=', num2str(A(s,1)), ', A2=', num2str(A(s,2)), ', A3=', num2str(A(s,3)), ...
           ', T1=', num2str(T(s,1)), ', T2=', num2str(T(s,2)), ', T3=', num2str(T(s,3))]);
    
    % Add legend
    legend({'Current (A)', 'Voltage (V)'}, 'Location', 'best');
    
    %% Construct W Matrix for Scenario s
    W = zeros(length(t), n);  % Initialize W matrix
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
    W_all{s} = W; % Store W matrix
end

%% Combine W and y Across All Scenarios for Joint DRT Estimation
W_combined = [];
y_combined = [];

for s = 1:num_scenarios
    W_combined = [W_combined; W_all{s}];
    y_adjusted = V_sd_all(s, :)' - OCV - R0 * ik_scenarios(s, :)';
    y_combined = [y_combined; y_adjusted];
end

%% Regularized Least Squares to Estimate Gamma
% Solve (W'W + lambda L'L) gamma = W'y
gamma_estimated = (W_combined' * W_combined + lambda * (L' * L)) \ (W_combined' * y_combined);
gamma_estimated(gamma_estimated < 0) = 0;  % Enforce non-negativity

% Normalize gamma_estimated for comparison
gamma_estimated_normalized = gamma_estimated / max(gamma_estimated);

%% Plot the DRT Comparison
figure;
plot(theta_discrete, gamma_discrete_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True DRT');
hold on;
plot(theta_discrete, gamma_estimated_normalized, 'b--', 'LineWidth', 2, 'DisplayName', 'Estimated DRT');
xlabel('ln(\tau) = \theta');
ylabel('\gamma(\theta)');
title(['Comparison of True DRT and Estimated DRT (\lambda = ', num2str(lambda), ')']);
legend('Location', 'Best');
grid on;
hold off;

% Adjust graph layout
set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);  % Screen ratio

%% Display Total Squared Voltage Error for Each Scenario
total_test_error = 0;
for s = 1:num_scenarios
    % Predict V_pred for the scenario
    V_pred = W_all{s} * gamma_estimated + R0 * ik_scenarios(s, :)' + OCV;
    
    % Actual noisy voltage
    V_actual = V_sd_all(s, :)';
    
    % Compute squared error
    error = sum((V_actual - V_pred).^2);
    total_test_error = total_test_error + error;
    
    fprintf('Scenario %d: Squared Voltage Error = %.4f\n', s, error);
end
fprintf('Total Squared Voltage Error across all scenarios: %.4f\n', total_test_error);
