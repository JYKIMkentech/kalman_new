clc; clear; close all;

% Set the random seed to ensure reproducibility
rng(0);  % You can change this seed to any number to generate a different fixed random sequence

% Parameters
num_scenarios = 10;  % Number of current scenarios
num_waves = 3;       % Number of sine waves per scenario
t = linspace(0, 1000, 10000);       % Time vector, 더 긴 시간 범위로 설정

% Constraints
T_min = 15;           % Minimum period value (approximate τ_min * 2π)
T_max = 250;          % Maximum period value (approximate τ_max * 2π)

% Initialize matrices for amplitudes, periods, and current scenarios
A = zeros(num_scenarios, num_waves);   % Amplitudes matrix
T = zeros(num_scenarios, num_waves);   % Periods matrix
ik_scenarios = zeros(num_scenarios, length(t)); % Current scenarios matrix

% Generate random amplitudes, periods, and current scenarios
for s = 1:num_scenarios
    % Random amplitudes that sum to 3
    temp_A = rand(1, num_waves);       % Generate 3 random amplitudes
    A(s, :) = 3 * temp_A / sum(temp_A);  % Normalize to make sum equal to 3
    
    % Random periods between T_min and T_max on a logarithmic scale
    log_T_min = log10(T_min);
    log_T_max = log10(T_max);
    T_log = log_T_min + (log_T_max - log_T_min) * rand(1, num_waves);
    T(s, :) = 10.^T_log;
    
    % Generate the current scenario as the sum of three sine waves
    ik_scenarios(s, :) = A(s,1)*sin(2*pi*t / T(s,1)) + ...
                         A(s,2)*sin(2*pi*t / T(s,2)) + ...
                         A(s,3)*sin(2*pi*t / T(s,3));
end

% Save the generated amplitudes, periods, and current scenarios to AS1.mat
save('AS1.mat', 'A', 'T', 'ik_scenarios', 't');

% Display the generated values for verification
disp('Amplitudes (A):');
disp(A);
disp('Periods (T):');
disp(T);

% Plot the 10 current scenarios in a 5x2 subplot grid
figure;
for s = 1:num_scenarios
    subplot(5, 2, s);  % Create a 5x2 grid of subplots
    plot(t, ik_scenarios(s, :), 'LineWidth', 1.5);
    title(['Scenario ', num2str(s)]);
    xlabel('Time (s)');
    ylabel('Current (A)');
    grid on;
end

% Adjust the layout for better spacing between subplots
sgtitle('Current Scenarios for 10 Randomized Cases');


