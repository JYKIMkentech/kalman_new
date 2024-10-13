clc; clear; close all;

%% 1. Parameters
n = 21;  % Number of RC elements
t = 0:0.01:100;  % Time vector (seconds)
dt = t(2) - t(1);  % Time step size
num_scenarios = 10;  % Number of current scenarios
lambda = 0.0409;  % Regularization parameter

%% 2. Define Amplitudes and Periods for Current Synthesis
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

%% 3. Generate Synthetic Current Data (Multi-Sine Approach)
ik_scenarios = zeros(num_scenarios, length(t)); % Initialize current scenarios

for s = 1:num_scenarios
    % Sum of three sine waves for each scenario
    ik_scenarios(s, :) = A(s,1)*sin(2*pi*t / T(s,1)) + ...
                         A(s,2)*sin(2*pi*t / T(s,2)) + ...
                         A(s,3)*sin(2*pi*t / T(s,3));
end

%% 4. True DRT Parameters (R_discrete)
mu = 10;
sigma = 5;
tau_discrete_linear = linspace(0.01, 20, n);  % Original linear tau values

% Normalize R_discrete_true to have a maximum value of 1
R_discrete_true_linear = normpdf(tau_discrete_linear, mu, sigma);
R_discrete_true_linear = R_discrete_true_linear / max(R_discrete_true_linear);

%% 5. DRT Using Logarithmic Time Scale (Natural Log)
% Define tau_j on a logarithmic scale using natural logarithms
a = 0.01;        % Minimum tau
b = 20;          % Maximum tau
theta_j = linspace(log(a), log(b), n); % Linearly spaced in log domain
tau_discrete_log = exp(theta_j);         % Exponentiate to get tau_j

% True DRT for log scale (assuming the same R_discrete_true)
R_discrete_true_log = normpdf(tau_discrete_log, mu, sigma);
R_discrete_true_log = R_discrete_true_log / max(R_discrete_true_log);  % Normalize

% Define Regularization Matrix L (1st order difference)
L = zeros(n-1, n);
for i = 1:n-1
    L(i, i) = -1;
    L(i, i+1) = 1;
end

%% 6. DRT Estimation for Each Scenario Using Logarithmic Time Scale
R_estimated_analytic = zeros(num_scenarios, n);   % Analytical estimates
R_estimated_quadprog = zeros(num_scenarios, n);  % Quadprog estimates
R_estimated_fmincon = zeros(num_scenarios, n);   % Fmincon estimates

% Initialize storage for V_est and V_sd
V_est_all = zeros(num_scenarios, length(t));  % Estimated voltages
V_sd_all = zeros(num_scenarios, length(t));   % Synthetic measured voltages

%% 7. Voltage Synthesis and Plotting Setup
figure(1);  
sgtitle('Current and Voltage for Each Scenario');

% Plot Current for each scenario
for s = 1:num_scenarios
    subplot(5, 2, s);
    yyaxis left
    plot(t, ik_scenarios(s, :), 'b-', 'LineWidth', 1.5);
    ylabel('Current (A)');
    xlabel('Time (s)');
    grid on;
end

%% 8. Processing Each Scenario
for s = 1:num_scenarios
    fprintf('Processing Scenario %d/%d...\n', s, num_scenarios);
    
    % Current for the scenario
    ik = ik_scenarios(s, :);  % Current scenario input
    
    %% 8.1. Initialize Voltage
    V_est = zeros(1, length(t));  % Model voltage calculated via n-RC model
    R0 = 0.1;  % Ohmic resistance (Ohms)
    OCV = 3.7; % Open Circuit Voltage (V)
    V_RC = zeros(n, length(t));  % RC voltages for each element
    
    %% 8.2. Initial Voltage Calculation (first time step)
    for i = 1:n
        V_RC(i, 1) = ik(1) * R_discrete_true_log(i) * (1 - exp(-dt / tau_discrete_log(i)));
    end
    V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));
    
    %% 8.3. Voltage Calculation for t > 1
    for k_idx = 2:length(t)
        for i = 1:n
            % Calculate RC voltage using previous time step data
            V_RC(i, k_idx) = exp(-dt / tau_discrete_log(i)) * V_RC(i, k_idx-1) + ...
                              R_discrete_true_log(i) * (1 - exp(-dt / tau_discrete_log(i))) * ik(k_idx);       
        end
        V_est(k_idx) = OCV + R0 * ik(k_idx) + sum(V_RC(:, k_idx));
    end
    
    % Store V_est for the current scenario
    V_est_all(s, :) = V_est;  % Save the calculated V_est for this scenario
    
    %% 8.4. Add Noise to the Voltage
    rng(s);  % Ensure reproducibility of noise for each scenario
    noise_level = 0.01;
    V_sd = V_est + noise_level * randn(size(V_est));  % V_sd = synthetic measured voltage
    
    % Store V_sd for the current scenario
    V_sd_all(s, :) = V_sd;  % Save the noisy V_sd for this scenario
    
    %% 8.5. Construct W Matrix for Log Scale
    W = zeros(length(t), n);  % Initialize W matrix
    V_RC = zeros(n, length(t));  % Reset V_RC for calculation
    
    for j = 1:n
        for k_idx = 1:length(t)
            if k_idx == 1
                W(k_idx, j) = ik(k_idx) * (1 - exp(-dt / tau_discrete_log(j)));
            else
                W(k_idx, j) = exp(-dt / tau_discrete_log(j)) * W(k_idx-1, j) + ...
                              ik(k_idx) * (1 - exp(-dt / tau_discrete_log(j)));
            end
        end
    end
    
    %% 8.6. Analytical Solution with Regularization
    % Compute (W^T W + lambda * L^T L)
    A_matrix = W' * W + lambda * (L' * L);
    
    % Compute W^T (V_sd - OCV - I * R0)
    b_vector = W' * (V_sd - OCV - ik * R0)';
    
    % Solve for R using the analytical solution
    R_analytic = A_matrix \ b_vector;
    
    % Enforce non-negativity
    R_analytic(R_analytic < 0) = 0;
    
    % Store Analytical DRT
    R_estimated_analytic(s, :) = R_analytic';
    
    %% 8.7. Quadratic Programming Solution
    H_qp = 2 * (W' * W + lambda * (L' * L));
    f_qp = -2 * W' * (V_sd - OCV - ik * R0)';
    
    % No equality constraints, only inequality (R >= 0)
    A_qp = [];
    b_qp = [];
    Aeq_qp = [];
    beq_qp = [];
    lb_qp = zeros(n,1); % R >= 0
    ub_qp = [];          % No upper bound
    
    % Options for quadprog
    options_qp = optimoptions('quadprog','Display','off');
    
    % Solve using quadprog
    R_quadprog = quadprog(H_qp, f_qp, A_qp, b_qp, Aeq_qp, beq_qp, lb_qp, ub_qp, [], options_qp);
    
    % Store Quadprog DRT
    R_estimated_quadprog(s, :) = R_quadprog';
    
    %% 8.8. Nonlinear Optimization using fmincon
    % Define the objective function
    objective = @(R) norm(V_sd - (OCV + ik * R0 + W * R), 2)^2 + lambda * norm(L * R', 2)^2;
    
    % Initial guess for R
    R_initial = zeros(n,1);
    
    % Define lower bounds (R >= 0)
    lb_fmin = zeros(n,1);
    ub_fmin = [];
    
    % Options for fmincon
    options_fmin = optimoptions('fmincon','Display','off','Algorithm','interior-point');
    
    % Solve using fmincon
    R_fmincon = fmincon(objective, R_initial, [], [], [], [], lb_fmin, ub_fmin, [], options_fmin);
    
    % Enforce non-negativity
    R_fmincon(R_fmincon < 0) = 0;
    
    % Store fmincon DRT
    R_estimated_fmincon(s, :) = R_fmincon';
    
    %% 8.9. Plot Voltage on Existing Subplots
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
end

%% 9. Plot the DRT Comparison for Each Scenario
for s = 1:num_scenarios
    figure(1 + s);  % DRT Comparison Figure for each scenario
    hold on;
    
    % Plot True DRT
    plot(tau_discrete_log, R_discrete_true_log, 'k-', 'LineWidth', 2, 'DisplayName', 'True DRT');
    
    % Plot Analytical DRT
    plot(tau_discrete_log, R_estimated_analytic(s, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical DRT');
    
    % Plot Quadratic Programming DRT
    plot(tau_discrete_log, R_estimated_quadprog(s, :), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Quadprog DRT');
    
    % Plot Nonlinear Optimization DRT
    plot(tau_discrete_log, R_estimated_fmincon(s, :), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Fmincon DRT');
    
    hold off;
    xlabel('\tau (Time Constant)');
    ylabel('R (Resistance)');
    title(['DRT Comparison for Scenario ', num2str(s), ' (\lambda = ', num2str(lambda), ')']);
    legend('Location', 'BestOutside');
    grid on;
end

%% 10. Summary Plot Comparing All Optimization Methods Across All Scenarios
figure;
hold on;
colors = lines(num_scenarios);

for s = 1:num_scenarios
    plot(tau_discrete_log, R_estimated_analytic(s, :), '--', 'Color', colors(s,:), 'LineWidth', 1);
    plot(tau_discrete_log, R_estimated_quadprog(s, :), ':', 'Color', colors(s,:), 'LineWidth', 1);
    plot(tau_discrete_log, R_estimated_fmincon(s, :), '-.', 'Color', colors(s,:), 'LineWidth', 1);
end

% Plot True DRT
plot(tau_discrete_log, R_discrete_true_log, 'k-', 'LineWidth', 2, 'DisplayName', 'True DRT');

xlabel('\tau (Time Constant)');
ylabel('R (Resistance)');
title('DRT Estimation Comparison Across All Scenarios');
legend({'Analytical DRT', 'Quadprog DRT', 'Fmincon DRT', 'True DRT'}, 'Location', 'BestOutside');
grid on;
hold off;

%% 11. Verify Dimensions
disp('Dimensions of matrices:');
disp(['W: ', mat2str(size(W))]);
disp(['R_analytic: ', mat2str(size(R_estimated_analytic))]);
disp(['R_quadprog: ', mat2str(size(R_estimated_quadprog))]);
disp(['R_fmincon: ', mat2str(size(R_estimated_fmincon))]);
disp(['L: ', mat2str(size(L))]);
