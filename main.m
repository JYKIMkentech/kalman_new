clear;
clc;

% Configuration parameters
Config.dt = 1;
Config.ik = 0.41; % Discharge current [A]
Config.R0 = 0.001884314;
Config.R1 = 0.045801322;
Config.C1 = 4846.080679;
Config.cap = (4.32019 * 3600) / 100; % nominal capacity [Ah]
Config.coulomb_efficient = 1;

% Initial conditions
SOC_initial = 1;
V1_initial = Config.ik * Config.R1 * (1 - exp(-Config.dt / (Config.R1 * Config.C1)));

% EKF initialization
P = eye(2); % Initial estimation error covariance

% Preallocate arrays for SOC, V1, and Vt
num_steps = 10; % Define the number of time steps for the simulation
intVar.SOC = zeros(1, num_steps);
intVar.V1 = zeros(1, num_steps);
intVar.Vt = zeros(1, num_steps);

% Initialize with the initial state
intVar.SOC(1) = SOC_initial;
intVar.V1(1) = V1_initial;
intVar.Vt(1) = ocv_soc(SOC_initial) - Config.R1 * V1_initial - Config.R0 * Config.ik;

% Simulation loop
for k = 2:num_steps
    % True SOC and Vt calculation
    [intVar.SOC(k), intVar.Vt(k)] = battery_model(intVar.SOC(k-1), intVar.V1(k-1), Config.ik, Config);
    
    % SOC estimation
    [intVar.SOC(k), intVar.V1(k), intVar.Vt(k), P] = soc_estimation(intVar.SOC(k-1), intVar.V1(k-1), intVar.Vt(k), Config.ik, Config, P);
end

% Plot SOC
figure;
plot(1:num_steps, intVar.SOC, 'b', 'LineWidth', 1.5); hold on;
plot(1:num_steps, intVar.SOC, 'r--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('SOC');
title('True SOC vs Estimated SOC');
legend('True SOC', 'Estimated SOC');
grid on;

% Plot Vt
figure;
plot(1:num_steps, intVar.Vt, 'b', 'LineWidth', 1.5); hold on;
plot(1:num_steps, intVar.Vt, 'r--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('V_t (V)');
title('True V_t vs Estimated V_t');
legend('True V_t', 'Estimated V_t');
grid on;
