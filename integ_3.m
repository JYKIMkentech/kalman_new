% 초기화
clear; clc; close all;

% 데이터 로드
load('C:\Users\deu04\OneDrive\바탕 화면\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
udds_current = meas.Current;
udds_voltage = meas.Voltage;
udds_time = meas.Time;

% SOC-OCV 로드
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);
ocv_values = soc_ocv(:, 2);

% 초기 추정값
initial_soc = 0.9901;
initial_params = [0.0254, 0.0147,60.474645]; % 초기 R0, R1, C 값

% 설정
Config.dt = mean(diff(udds_time));
Config.cap = 2.90; % Ah
Config.coulomb_efficient = 1;

% 비용 함수 정의
cost_function = @(params) func_cost(params, initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values);

% 최적화 옵션 설정
options = optimoptions('fmincon', 'Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-8); 

% 제약 조건 설정 (비음수 제약 조건 추가)
lb = [0, 0, 0]; % 하한 (비음수)
ub = []; % 상한 (없음)

% 최적화 실행
[optimized_params, resnorm] = fmincon(cost_function, initial_params, [], [], [], [], lb, ub, [], options);

% 결과 출력
fprintf('Optimized R0: %.6f\n', optimized_params(1));
fprintf('Optimized R1: %.6f\n', optimized_params(2));
fprintf('Optimized C : %.6f\n', optimized_params(3));

% 최적화된 값으로 결과 플롯
[~, Vt_est] = voltage_error(optimized_params, initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values);

figure;
plot(udds_time, udds_voltage, 'b', 'LineWidth', 1.5); hold on;
plot(udds_time, Vt_est, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Measured vs Estimated Terminal Voltage during UDDS Cycle');
legend('Measured V_t', 'Estimated V_t');
grid on;

% Cost function
function cost = func_cost(params, initial_soc, current, voltage, Config, soc_values, ocv_values)
    residuals = voltage_error(params, initial_soc, current, voltage, Config, soc_values, ocv_values);
    cost = sum(residuals.^2);
end

% Vt-Vest residuals 
function [residuals, Vt_est] = voltage_error(params, initial_soc, current, voltage, Config, soc_values, ocv_values)
    R0 = params(1);
    R1 = params(2);
    C1 = params(3);

    num_samples = length(current);
    SOC_est = zeros(num_samples, 1);
    V1_est = zeros(num_samples, 1);
    Vt_est = zeros(num_samples, 1);

    SOC_est(1) = initial_soc;
    V1_est(1) = current(1) * R1 * (1 - exp(-Config.dt / (R1 * C1)));
    Vt_est(1) = voltage(1);

    for k = 2:num_samples
        ik = current(k);
        SOC_est(k) = SOC_est(k-1) + (Config.dt / (Config.cap * 3600)) * Config.coulomb_efficient * ik;
        V1_est(k) = exp(-Config.dt / (R1 * C1)) * V1_est(k-1) + (1 - exp(-Config.dt / (R1 * C1))) * ik * R1;
        Vt_est(k) = interp1(soc_values, ocv_values, SOC_est(k), 'linear', 'extrap') + V1_est(k) + R0 * ik;
    end

    residuals = Vt_est - voltage;
end
