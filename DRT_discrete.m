clc; clear; close all;

% Parameters
n = 5;  % RC 요소의 개수
m = 21;  % 이산적인 tau와 R 값의 개수 (n과 다름)
t = 0:0.01:100;  % 시간 벡터 (이산 시간)
dt = t(2) - t(1);

% 합성 전류 데이터 (사인파들의 합)
A = 1; 
T1 = 1;
T2 = 5;
T3 = 20;
I1 = A * sin(2 * pi * t / T1);
I2 = A * sin(2 * pi * t / T2);
I3 = A * sin(2 * pi * t / T3);
ik = I1 + I2 + I3;  % 총 전류

% tau에 대한 정규분포 파라미터
mu = 10;
sigma = 5;
tau_discrete = linspace(0.01, 20, m);  % 이산적인 tau 값들 (m 값)

% tau_discrete에 해당하는 R_discrete을 normpdf로 계산
R_discrete = normpdf(tau_discrete, mu, sigma);

% R_discrete을 정규화하여 최대 값이 1이 되도록 조정
R_discrete = R_discrete / max(R_discrete);  % R_discrete의 최대 값은 이제 1

% Initialize voltage
V_est = zeros(1, length(t));  % Estimated voltage
R0 = 0.1;  % Internal resistance
OCV = 0;   % Open circuit voltage
V_RC = zeros(n, length(t));  % RC voltages for each element

%% Initial voltage calculation (first time step)
for i = 1:n
    V_RC(i, 1) = ik(1) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i))); % ik(1) = 0 
    V_est(1) = OCV + R0 * ik(1) + V_RC(i, 1);
end

%% Discrete-time voltage calculation (n-RC system for subsequent time steps)
for k = 2:length(t)
    for i = 1:n
        % Calculate RC voltages based on previous time step
        V_RC(i, k) = exp(-dt / tau_discrete(i)) * V_RC(i, k-1) + R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k);       
    end
    V_est(k) = OCV + R0 * ik(k) + sum(V_RC(:, k));
end

%% Add noise to the voltage
noise_level = 0.01;
V_noisy = V_est + noise_level * randn(size(V_est));

%% Optimization using fmincon to find the best R values
initial_R = R_discrete;% ones(1, m)/2;  %R_discrete   % Initial guess for R_discrete
lb = [];%zeros(1, m);  % Lower bound (R cannot be negative)
ub = [];%ones(1, m);   % Upper bound (normalize to 1)

% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
R_optimized = fmincon(@cost_function, initial_R, [], [], [], [], lb, ub, [], options, tau_discrete, ik, V_noisy, dt, n, R0, OCV, t);

%% Plot the DRT comparison (True DRT vs Optimized DRT)
figure;
plot(tau_discrete, R_discrete, 'b-', 'LineWidth', 2);  % True DRT (blue solid line)
hold on;
stem(tau_discrete, R_discrete, 'bo', 'LineWidth', 1.5);  % True DRT points (blue circles)

plot(tau_discrete, R_optimized, 'r-', 'LineWidth', 2);  % Optimized DRT (red solid line)
stem(tau_discrete, R_optimized, 'ro', 'LineWidth', 1.5);  % Optimized DRT points (red circles)

xlabel('\tau (Time Constant)');
ylabel('R (Resistance)');
legend('True DRT', 'Optimized DRT');
title('True vs Optimized Distribution of Relaxation Times (DRT)');
grid on;

%% Functions

% Residuals function for calculating V_est with given R_discrete
function V_est = calculate_voltage(R_discrete, tau_discrete, ik, dt, n, R0, OCV, t)
    V_est = zeros(1, length(t));  % Initialize estimated voltage
    V_RC = zeros(n, length(t));  % RC voltages for each element

    % Initial voltage calculation (first time step)
    for i = 1:n
        V_RC(i, 1) = ik(1) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i)));
        V_est(1) = OCV + R0 * ik(1) + V_RC(i, 1);
    end

    % Discrete-time voltage calculation for subsequent time steps
    for k = 2:length(t)
        for i = 1:n
            V_RC(i, k) = exp(-dt / tau_discrete(i)) * V_RC(i, k-1) + R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k);
        end
        V_est(k) = OCV + R0 * ik(k) + sum(V_RC(:, k));
    end
end

% Cost function for fmincon (sum of squared residuals)
function cost = cost_function(R_discrete, tau_discrete, ik, V_noisy, dt, n, R0, OCV, t)
    % Calculate the estimated voltage for the current R_discrete
    V_est = calculate_voltage(R_discrete, tau_discrete, ik, dt, n, R0, OCV, t);
    
    % Compute the sum of squared differences (residuals)
    residuals = V_noisy - V_est;
    cost = sum(residuals.^2);
end



