clc; clear; close all;

%% Parameters
n = 21;  % Number of RC elements (extendable to any order)
t = 0:0.01:100;  % Time vector (discrete time)

% Synthetic current data (sum of sine waves)
A = 1; 
T1 = 1;
T2 = 5;
T3 = 20;
I1 = A * sin(2 * pi * t / T1);
I2 = A * sin(2 * pi * t / T2);
I3 = A * sin(2 * pi * t / T3);
ik = I1 + I2 + I3;  % Total current

% Parameters for the normal distribution of tau
mu = 10;
sigma = 5;
tau_discrete = linspace(0.01, 20, n);  % Discrete tau values for RC elements

% Use normpdf to directly calculate the R_discrete corresponding to tau_discrete
R_discrete = normpdf(tau_discrete, mu, sigma);

% Normalize R_discrete so that the maximum value is 1
R_discrete = R_discrete / max(R_discrete);  % Now the maximum value of R_discrete is 1

% Initialize voltage
V_est = zeros(1, length(t));  % Estimated voltage
R0 = 0.1;  % Internal resistance
OCV = 0;   % Open circuit voltage
V_RC = zeros(n, length(t));  % RC voltages for each element

%% Voltage calculation for each time step
for k = 1:length(t)
    if k < length(t)
        dt = t(k+1) - t(k);  % Calculate dynamic time step for each time
    else
        dt = t(k) - t(k-1);  % For the last step, use the previous dt
    end
    
    % Compute RC voltages for each RC element
    for i = 1:n
        if k == 1
            % First time step, initial V_RC calculation
            V_RC(i, k) = ik(k) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i)));
        else
            % Subsequent time steps, accumulate V_RC
            V_RC(i, k) = exp(-dt / tau_discrete(i)) * V_RC(i, k-1) + ...
                         R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k);
        end
    end

    % Final voltage summation: OCV + R0 * I + sum(V_RC)
    V_est(k) = OCV + R0 * ik(k) + sum(V_RC(:, k));
end

%% Add noise to the voltage
noise_level = 0.01;
V_noisy = V_est + noise_level * randn(size(V_est));

%% Plot the results
figure;

% Subplot for current
subplot(3,1,1);
plot(t, ik, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current (A)');
title('Current');
grid on;

% Subplot for voltage
subplot(2,1,2);
plot(t, V_est, 'b-', 'LineWidth', 1.5); hold on;
plot(t, V_noisy, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_{est}', 'V_{noisy}');
title('Estimated Voltage and Noisy Voltage (n-RC System)');
grid on;

% Subplot for tau vs R (True DRT)
figure(2)

plot(tau_discrete, R_discrete, 'b-', 'LineWidth', 1.5);  % Smooth curve
hold on;
stem(tau_discrete, R_discrete, 'r', 'LineWidth', 1.5);  % Vertical lines
plot(tau_discrete, R_discrete, 'ro', 'LineWidth', 1.5);  % Points on the curve
xlabel('\tau (Time Constant)');
ylabel('R (Resistance)');
title('Distribution of Relaxation Times (DRT)');
grid on;

