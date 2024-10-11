clc; clear; close all;

%% parameters 

n = 21;  % Number of RC elements
t = 0:0.01:100;  % Time vector
dt = t(2) - t(1);
num_scenarios = 10;  % Number of current scenarios
lambda = 0.1;  % Regularization parameter

% synthetic current parameters 

Amp = linspace(1, 10, num_scenarios);  % Amplitude
T = [1, 2, 5, 10, 20, 25, 30, 35, 40, 50];  % Period
ik_scenarios = zeros(num_scenarios, length(t));  

% True DRT parameters (R_discrete)
mu = 10;
sigma = 5;
tau_discrete = linspace(0.01, 20, n);  % Discrete tau values

% 1차 미분 L 행렬
L = zeros(n-1, n);
for i = 1:n-1
    L(i, i) = -1;
    L(i, i+1) = 1;
end

%% DRT 

% True DRT [Tau,R]
R_discrete_true = normpdf(tau_discrete, mu, sigma);
R_discrete_true = R_discrete_true / max(R_discrete_true);  % Normalize to max value of 1

% DRT 저항 초기값값 [Tau,R]
R_analytical_all = zeros(num_scenarios, n);  % Analytical DRT estimates

%% 전류 합성

% 10개 전류 시나리오
for k = 1:num_scenarios
    ik_scenarios(k, :) = Amp(k) * sin(2 * pi * t / T(k));
end


%% 전압 합성
figure(1);  
% 매 시나리오마다 전류,전압 합성 반복
for s = 1:num_scenarios
    fprintf('Processing Scenario %d/%d...\n', s, num_scenarios);
    
    % Current for the scenario
    ik = ik_scenarios(s, :); % 전류 시나리오 입력
    
    %% Initialize voltage
    V_est = zeros(1, length(t));  % n-RC model을 통해 계산된 model voltage
    R0 = 0.1;  % Ohmic resistance
    OCV = 0;   % OCV
    V_RC = zeros(n, length(t));  % RC voltages for each element
    
    %% 초기 전압 계산  (first time step)
    for i = 1:n
        V_RC(i, 1) = ik(1) * R_discrete_true(i) * (1 - exp(-dt / tau_discrete(i)));
    end
    V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));
    
    %% time > 1 이후부터 전압 계산 
    % 계산된 V_est = OCV  + IR0 + V_RC

    for k_idx = 2:length(t)
        for i = 1:n
            % Calculate RC voltages based on previous time step
            V_RC(i, k_idx) = exp(-dt / tau_discrete(i)) * V_RC(i, k_idx-1) + ...
                             R_discrete_true(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k_idx);       
        end
        V_est(k_idx) = OCV + R0 * ik(k_idx) + sum(V_RC(:, k_idx));
    end
    
    %% Add noise to the voltage
    rng(0);  % 노이즈 일정하게 추가
    noise_level = 0.01; %
    V_sd = V_est + noise_level * randn(size(V_est)); % Vsd = synthethic data voltage = 측정 전압 (such as UDDS) 
    
    %% COST = (측정 전압 - 모델 전압 ) ^2 + penalty 
    % cost f(R) = (Vsd - Vest)^2 + lambda * ( LR ) ^2 = (Vsd - W * R )^2 + lambda * ( LR ) ^2 
    
    W = zeros(length(t), n);  % Initialize W matrix
    for k_idx = 1:length(t)
        for i = 1:n
            if k_idx == 1 % 1초일때 행렬 
                W(k_idx, i) = ik(k_idx) * (1 - exp(-dt / tau_discrete(i)));
            else % 1초 이상일때 W 행렬
                W(k_idx, i) = exp(-dt / tau_discrete(i)) * W(k_idx-1, i) + ...
                              ik(k_idx) * (1 - exp(-dt / tau_discrete(i)));
            end
        end
    end
    
    %% 정규화을 통한 Analytical solution 

    R_analytical = (W' * W + lambda * (L' * L)) \ (W' * V_sd'); 
    R_analytical(R_analytical < 0) = 0;  % Enforce non-negativity
    
    %% 결과값 저장
    R_analytical_all(s, :) = R_analytical';
    
    %% 전류 전압 subplot 그래프
    subplot(5, 2, s);
    yyaxis left
    plot(t, ik, 'b-', 'LineWidth', 1.5);
    ylabel('Current (A)');
    ylim([min(ik)-1, max(ik)+1]);
    
    yyaxis right
    plot(t, V_sd, 'r-', 'LineWidth', 1.5);
    ylabel('Voltage (V)');
    ylim([min(V_sd)-0.1, max(V_sd)+0.1]);
    
    title(['Scenario ', num2str(s), ': A=', num2str(Amp(s)), ', T=', num2str(T(s))]);
    xlabel('Time (s)');
    grid on;
end

%% Enhance the Current and Voltage Figure
figure(1);
sgtitle('Current and Voltage for Each Scenario');

%% Plot the DRT comparison for each scenario in separate figures
for s = 1:num_scenarios
    figure(2 + s);  % DRT Comparison Figure for each scenario
    hold on;
    
    % Plot True DRT
    plot(tau_discrete, R_discrete_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True DRT');
    
    % Plot Analytical DRT
    plot(tau_discrete, R_analytical_all(s, :), ':', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'Analytical DRT');
    
    hold off;
    xlabel('\tau (Time Constant)');
    ylabel('R (Resistance)');
    title(['DRT Comparison for Scenario ', num2str(s), ' (\lambda = ', num2str(lambda), ')']);
    legend('True DRT', 'Analytical DRT', 'Location', 'BestOutside');
    grid on;
end

%% Functions

% % Residuals function for calculating V_est with given R_discrete
% function V_est = calculate_voltage(R_discrete, tau_discrete, ik, dt, n, R0, OCV, t)
%     V_est = zeros(1, length(t));  % Initialize estimated voltage
%     V_RC = zeros(n, length(t));  % RC voltages for each element
% 
%     % Initial voltage calculation (first time step)
%     for i = 1:n
%         V_RC(i, 1) = ik(1) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i)));
%     end
%     V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));
% 
%     % Discrete-time voltage calculation for subsequent time steps
%     for k_idx = 2:length(t)
%         for i = 1:n
%             V_RC(i, k_idx) = exp(-dt / tau_discrete(i)) * V_RC(i, k_idx-1) + ...
%                              R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k_idx);
%         end
%         V_est(k_idx) = OCV + R0 * ik(k_idx) + sum(V_RC(:, k_idx));
%     end
% end










