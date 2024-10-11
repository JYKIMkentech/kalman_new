clc; clear; close all;

%% Parameters
n = 21;  % Number of RC elements
t = 0:0.01:100;  % Time vector
dt = t(2) - t(1);
num_scenarios = 10;  % Total number of current scenarios

% Construct the first derivative matrix L for regularization
L = zeros(n-1, n);
for i = 1:n-1
    L(i, i) = -1;
    L(i, i+1) = 1;
end

% Synthetic current data (sum of sine waves for multiple scenarios)
Amp = linspace(1, 10, num_scenarios);  % Amplitude values from 1 to 10
T = [1, 2, 5, 10, 20, 25, 30, 35, 40, 50];  % Different period values for each scenario
ik_scenarios = zeros(num_scenarios, length(t));  % Initialize current for each scenario

% Generate the 10 sine waves
for k = 1:num_scenarios
    ik_scenarios(k, :) = Amp(k) * sin(2 * pi * t / T(k));
end

% Parameters for the true DRT (R_discrete)
mu = 10;
sigma = 5;
tau_discrete = linspace(0.01, 20, n);  % Discrete tau values

% Calculate true R_discrete using a normal distribution
R_discrete_true = normpdf(tau_discrete, mu, sigma);
R_discrete_true = R_discrete_true / max(R_discrete_true);  % Normalize to max value of 1