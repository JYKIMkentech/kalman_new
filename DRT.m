clc; clear; close all;

%% Synthetic Current 
% Time vector
t = 0:0.01:100;  % Time from 0 to 100 seconds with a step of 0.01

% Synthetic current data (sum of sine waves)
A = 1; 
T1 = 1;        % Period for first sine wave
T2 = 5;        % Period for second sine wave
T3 = 20;       % Period for third sine wave

% Current components (synthetic current data)
I1 = A * sin(2 * pi * t / T1);
I2 = A * sin(2 * pi * t / T2);
I3 = A * sin(2 * pi * t / T3);
ik = I1 + I2 + I3;  % Total current

%% Normal PDF 가정 (DRT에 사용)
mu = 10;      % Mean of the distribution
sigma = 5;    % Standard deviation of the distribution

% Time constant (tau) range
tau = 0.1:0.1:20;  % tau starts from 0.1 instead of 0

% Normal PDF
pdf_tau = (1/(sigma * sqrt(2*pi))) * exp(-(tau - mu).^2 / (2 * sigma^2));
pdf_tau = pdf_tau / max(pdf_tau);  % Normalize the PDF so max value is 1

%% Discretize tau and corresponding R for 5-RC model
num_RC = 100;  % Number of RC circuits
tau_discrete = linspace(min(tau), max(tau), num_RC);  % Discrete tau values (starting > 0)

% R_discrete follows the normal PDF distribution
R_discrete = interp1(tau, pdf_tau, tau_discrete);  % Interpolate R values from the PDF

%% Calculate V_est using the discrete RC model (with initial R_discrete)
V_est = zeros(size(t));
R0 = 0.1;  % Example internal resistance
OCV = 0;   % Open circuit voltage

for i = 1:num_RC
    V_est = V_est + ik .* R_discrete(i) .* (1 - exp(-t / tau_discrete(i)));
end

V_est = OCV + ik .* R0 + V_est;
V_est(1) = 0;  % Set initial voltage

%% Add noise to the voltage
noise_level = 0.1;  % 1% noise
V_noisy = V_est + noise_level * randn(size(V_est));

%% Define cost and residual functions (fitting only R, not tau)

% Residual function: calculates the difference between V_noisy and V_est
residual_fn = @(R_params) V_noisy - calculate_V_est(R_params, tau_discrete, t, ik, R0, OCV);

% Cost function: sum of squared residuals with L2 regularization
lambda = 0.1;  % Regularization parameter
cost_fn = @(R_params) sum(residual_fn(R_params).^2) + lambda * sum(R_params.^2);

%% fmincon setup (fitting only R)
% Initial guess for R (fitting parameter)
initial_guess_R = ones(1, num_RC);

% Bounds for optimization (optional)
lb = zeros(1, num_RC);  % Lower bound for R
ub = ones(1, num_RC) * 10;  % Upper bound for R

% Call fmincon for optimization (fitting only R, tau is fixed)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1e-8, 'TolX', 1e-8);
[estimated_R, fval] = fmincon(cost_fn, initial_guess_R, [], [], [], [], lb, ub, [], options);

%% Plot the results with subplots (Current and Voltage)
figure;

% Subplot for current
subplot(2,1,1);  % 2 rows, 1 column, 1st subplot
plot(t, ik, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current (A)');
title('Current');
grid on;

% Subplot for voltage
subplot(2,1,2);  % 2 rows, 1 column, 2nd subplot
plot(t, V_noisy, 'r--', 'LineWidth', 1.5); hold on;
plot(t, calculate_V_est(estimated_R, tau_discrete, t, ik, R0, OCV), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_{noisy}', 'V_{estimated}');
title('Noisy Voltage and Optimized Estimated Voltage (Fitting R only)');
grid on;

%% Additional figure for (tau, R) comparison

figure;

% Subplot 1: Initial Normal PDF-based tau and R
subplot(2, 1, 1);
scatter(tau_discrete, R_discrete, 100, 'filled');  % Initial (tau, R) from normal PDF
xlabel('Time Constant (tau)');
ylabel('Resistance (R)');
title('Initial (tau, R) based on Normal PDF');
grid on;

% Subplot 2: Fitted R with fixed tau
subplot(2, 1, 2);
scatter(tau_discrete, estimated_R, 100, 'filled');  % Fitted R with fixed tau
xlabel('Time Constant (tau)');
ylabel('Resistance (R)');
title('Fitted (R) with Fixed (tau)');
grid on;

%% Helper function: Calculate V_est
function V_est = calculate_V_est(R_params, tau_vals, t, ik, R0, OCV)
    num_RC = length(R_params);
    V_est = OCV + ik .* R0;  % Start with OCV and R0 contribution
    
    % Add RC elements' contributions with fixed tau
    for i = 1:num_RC
        V_est = V_est + ik .* R_params(i) .* (1 - exp(-t / tau_vals(i)));
    end
end


