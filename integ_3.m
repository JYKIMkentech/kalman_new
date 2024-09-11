clear; clc; close all;

% 데이터 로드
load('C:\Users\USER\Desktop\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
udds_current = meas.Current;  % 전류 데이터
udds_voltage = meas.Voltage;  % 전압 데이터
udds_time = meas.Time;        % 시간 데이터

% SOC-OCV 데이터 로드
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);  % SOC 값
ocv_values = soc_ocv(:, 2);  % OCV 값

% 초기 추정값 설정
initial_soc = 0.9901;  % 초기 SOC 값
fixed_R0 = 0.025426;   % 고정된 R0 값
initial_params = [0.014184, 21.465086]; % 초기 R1, C 값

% 설정값
Config.dt = mean(diff(udds_time));  % 시간 간격
Config.cap = 2.90;                  % 배터리 용량 (Ah)
Config.coulomb_efficient = 1;       % 쿨롱 효율

% 비용 함수 정의 (시간 가중치 미적용)
cost_function = @(params) func_cost(params, fixed_R0, initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values);

% 최적화 옵션 설정
options = optimoptions('fmincon', 'Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-8); 

% 제약 조건 설정 (비음수 제약 조건 추가)
lb = [0, 0]; % 하한 (비음수)
ub = [];     % 상한 (없음)

% 최적화 실행
[optimized_params, resnorm] = fmincon(cost_function, initial_params, [], [], [], [], lb, ub, [], options);

% 최적화된 값 출력
fprintf('Optimized R1: %.6f\n', optimized_params(1));
fprintf('Optimized C : %.6f\n', optimized_params(2));

% 최적화된 값으로 결과 플롯
[~, Vt_est] = voltage_error([fixed_R0, optimized_params], initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values);

figure;
plot(udds_time, udds_voltage, 'b', 'LineWidth', 1.5); hold on;
plot(udds_time, Vt_est, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Measured vs Estimated Terminal Voltage during UDDS Cycle');
legend('Measured V_t', 'Estimated V_t');
xlim([0 100])
grid on;

% R1, C1 범위 설정
R1_range = linspace(0.001, 0.1, 1000);
C1_range = linspace(1, 100, 1000);

% 비용 함수 계산을 위한 배열 초기화
cost_values = zeros(length(R1_range), length(C1_range));

% 비용 함수 계산
for i = 1:length(R1_range)
    for j = 1:length(C1_range)
        params = [R1_range(i), C1_range(j)];
        cost_values(i, j) = func_cost(params, fixed_R0, initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values);
    end
end

% 3D 그래프 플롯
[R1_mesh, C1_mesh] = meshgrid(R1_range, C1_range);

figure;
surf(R1_mesh, C1_mesh, cost_values');
xlabel('R1 (Ohms)');
ylabel('C1 (Farads)');
zlabel('Cost');
title('Cost Function Surface for R0, R1, C1');
colormap('viridis');
colorbar;
grid on;

% 비용 함수 정의 (시간 가중치 없음)
function cost = func_cost(params, R0, initial_soc, current, voltage, Config, soc_values, ocv_values)
    % 잔차(residual) 계산
    [residuals, ~] = voltage_error([R0, params], initial_soc, current, voltage, Config, soc_values, ocv_values);
    
    % 비용 계산
    cost = sum(residuals.^2);
end

% Vt-Vest 잔차 계산 함수
function [residuals, Vt_est] = voltage_error(params, initial_soc, current, voltage, Config, soc_values, ocv_values)
    R0 = params(1);  % 고정된 R0
    R1 = params(2);  % 최적화된 R1
    C1 = params(3);  % 최적화된 C1

    num_samples = length(current);  % 샘플 개수
    SOC_est = zeros(num_samples, 1);  % SOC 추정값
    V1_est = zeros(num_samples, 1);   % V1 추정값
    Vt_est = zeros(num_samples, 1);   % 단자 전압 추정값

    SOC_est(1) = initial_soc;  % 초기 SOC 설정
    V1_est(1) = current(1) * R1 * (1 - exp(-Config.dt / (R1 * C1)));  % 초기 V1 추정
    Vt_est(1) = voltage(1);  % 초기 전압 설정

    for k = 2:num_samples
        ik = current(k);
        SOC_est(k) = SOC_est(k-1) + (Config.dt / (Config.cap * 3600)) * Config.coulomb_efficient * ik;
        V1_est(k) = exp(-Config.dt / (R1 * C1)) * V1_est(k-1) + (1 - exp(-Config.dt / (R1 * C1))) * ik * R1;
        Vt_est(k) = interp1(soc_values, ocv_values, SOC_est(k), 'linear', 'extrap') + V1_est(k) + R0 * ik;
    end

    residuals = Vt_est - voltage;  % 잔차 계산
end

