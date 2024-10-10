clc; clear; close all;

% Parameters
n = 21;  % Number of RC elements
t = 0:0.01:100;  % Time vector
dt = t(2) - t(1);

% Synthetic current data (sum of sine waves)
Amp = linspace(1, 10, 10);  % Amplitude values from 1 to 10 (arithmetic progression)
T = [1, 2, 5, 10, 20, 25, 30, 35, 40, 50];  % Different period values for each signal
ik_scenarios = zeros(10, length(t));  % Initialize current for each scenario

% Generate the 10 sine waves
for k = 1:10
    ik_scenarios(k, :) = Amp(k) * sin(2 * pi * t / T(k));
end

% Parameters for the true DRT (R_discrete)
mu = 10;
sigma = 5;
tau_discrete = linspace(0.01, 20, n);  % Discrete tau values

% Calculate true R_discrete using a normal distribution
R_discrete = normpdf(tau_discrete, mu, sigma);
R_discrete = R_discrete / max(R_discrete);  % Normalize to max value of 1

% Initialize arrays for voltage estimations
V_est_scenarios = zeros(10, length(t));  % Estimated voltage for each scenario
R0 = 0.1;  % Internal resistance
OCV = 0;   % Open circuit voltage
V_RC = zeros(n, length(t));  % RC voltages for each element

%% Loop over each scenario for current and voltage calculations
for s = 1:10
    ik = ik_scenarios(s, :);  % Get the current for the scenario
    
    % Initial voltage calculation (first time step)
    for i = 1:n
        V_RC(i, 1) = ik(1) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i)));
    end
    V_est_scenarios(s, 1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));
    
    % Discrete-time voltage calculation for subsequent time steps
    for k = 2:length(t)
        for i = 1:n
            V_RC(i, k) = exp(-dt / tau_discrete(i)) * V_RC(i, k-1) + R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k);
        end
        V_est_scenarios(s, k) = OCV + R0 * ik(k) + sum(V_RC(:, k));
    end
end

%% Analytical solution using normal equations (keep this active)
W = zeros(length(t), n);
for i = 1:n
    for k = 1:length(t)
        if k == 1
            W(k, i) = ik_scenarios(1, k) * (1 - exp(-dt / tau_discrete(i)));
        else
            W(k, i) = exp(-dt / tau_discrete(i)) * W(k-1, i) + ik_scenarios(1, k) * (1 - exp(-dt / tau_discrete(i)));
        end
    end
end

y = V_est_scenarios(1, :)' - OCV - R0 * ik_scenarios(1, :)';

% Compute A and b for normal equations
lambda = 0.1;  % Regularization parameter
L = diag(ones(1, n-1), 1) - diag(ones(1, n), 0);  % Derivative matrix
A_matrix = W' * W + lambda * (L' * L);  % Rename to A_matrix
b = W' * y;

% Solve the normal equations
R_analy = A_matrix \ b;

% Enforce non-negativity
R_analy(R_analy < 0) = 0;

%% Plotting current and voltage for all scenarios
figure;
for s = 1:10
    subplot(5, 2, s);
    
    % Plot current (ik) and voltage (V_est) in the same subplot
    yyaxis left
    plot(t, ik_scenarios(s, :), 'b-', 'LineWidth', 1.5);
    ylabel('Current (A)');
    
    yyaxis right
    plot(t, V_est_scenarios(s, :), 'r-', 'LineWidth', 1.5);
    ylabel('Voltage (V)');
    
    % Correctly set the title for each scenario
    title(['Scenario ', num2str(s), ': A=', num2str(Amp(s)), ', T=', num2str(T(s))]);
    xlabel('Time (s)');
    grid on;
end

%% Functions

% Residuals function for calculating V_est with given R_discrete
function V_est = calculate_voltage(R_discrete, tau_discrete, ik, dt, n, R0, OCV, t)
    V_est = zeros(1, length(t));  % Initialize estimated voltage
    V_RC = zeros(n, length(t));  % RC voltages for each element

    % Initial voltage calculation (first time step)
    for i = 1:n
        V_RC(i, 1) = ik(1) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i)));
    end
    V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));

    % Discrete-time voltage calculation for subsequent time steps
    for k = 2:length(t)
        for i = 1:n
            V_RC(i, k) = exp(-dt / tau_discrete(i)) * V_RC(i, k-1) + R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k);
        end
        V_est(k) = OCV + R0 * ik(k) + sum(V_RC(:, k));
    end
end

% Cost function for fmincon (sum of squared residuals with regularization)
function cost = cost_function(R_discrete, tau_discrete, ik, V_noisy, dt, n, R0, OCV, t, lambda, L)
    % Calculate the estimated voltage for the current R_discrete
    V_est = calculate_voltage(R_discrete, tau_discrete, ik, dt, n, R0, OCV, t);
    
    % Compute the sum of squared differences (residuals)
    residuals = V_noisy - V_est;
    data_fidelity = sum(residuals.^2);
    
    % Regularization term (first derivative penalty)
    regularization = lambda * norm(L * R_discrete', 2)^2;
    
    % Total cost
    cost = data_fidelity + regularization;
end


