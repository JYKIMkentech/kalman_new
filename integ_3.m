clear; clc; close all;

% 데이터 로드
load('C:\Users\deu04\OneDrive\바탕 화면\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
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

% 설정값
Config.cap = 2.90;                  % 배터리 용량 (Ah)
Config.coulomb_efficient = 1;       % 쿨롱 효율

% R1과 C1의 범위 설정
R1_range = linspace(0.1, 0.2, 6);  % R1 값의 범위 (Ohms)
C1_range = linspace(10, 1300, 6);  % C1 값의 범위 (Farads)

% 비용 함수 계산을 위한 배열 초기화
cost_values = zeros(length(R1_range), length(C1_range));

% 비용 함수 계산
for i = 1:length(R1_range)
    for j = 1:length(C1_range)
        params = [R1_range(i), C1_range(j)];
        cost_values(i, j) = func_cost(params, fixed_R0, initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values, udds_time);
    end
end

% 3D 그래프 플롯
[R1_mesh, C1_mesh] = meshgrid(R1_range, C1_range);

figure;
surf(R1_mesh, log10(C1_mesh), log10(cost_values'), 'EdgeColor', 'none');  
xlabel('R1 ');
ylabel('log C1 ');
zlabel('Cost');
title('log Cost');
colorbar;
grid on;

% 2D 등고선 그래프 플롯
figure;
contourf(R1_mesh, log10(C1_mesh), log10(cost_values'), 20, 'LineColor', 'none');  
ylabel('logC1 ');
title('log Cost ');
colorbar;
grid on;

% 초기 추정값 (초기 guess)
initial_params = [0.014184, 1.465];  % R1, C1 초기값

% fmincon 옵션 설정
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
lb = [];
ub = []; %initial_params * 3;

% fmincon 문제 정의
problem = createOptimProblem('fmincon', ...
                             'objective', @(params) func_cost(params, fixed_R0, initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values, udds_time), ...
                             'x0', initial_params, ...
                             'lb', lb, ...
                             'ub', ub, ...
                             'options', options);

ms = MultiStart('UseParallel', true, 'StartPointsToRun', 'all');

% MultiStart 실행
num_start_points = 20;  % 시도할 시작점의 수
[optimal_params, fval] = run(ms, problem, num_start_points);

% 최적화 결과
optimal_R1 = optimal_params(1);
optimal_C1 = optimal_params(2);

fprintf('Optimal R1: %.6f Ohms\n', optimal_R1);
fprintf('Optimal C1: %.6f Farads\n', optimal_C1);

% 최적화된 파라미터로 Vt_est 계산
[residuals, Vt_est] = voltage_error([fixed_R0, optimal_R1, optimal_C1], initial_soc, udds_current, udds_voltage, Config, soc_values, ocv_values, udds_time);

% 최적화 결과 시각화
figure;
plot(udds_time, udds_voltage, 'b', 'LineWidth', 1.5); hold on;
plot(udds_time, Vt_est, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Measured vs Estimated Terminal Voltage during UDDS Cycle');
legend('Measured V_t', 'Estimated V_t');
xlim([0 100]);  % 필요에 따라 범위 조정
grid on;

% 비용 함수 정의 (시간 가중치 없음)
function cost = func_cost(params, R0, initial_soc, current, voltage, Config, soc_values, ocv_values, time_data)
    % 잔차(residual) 계산
    [residuals, ~] = voltage_error([R0, params], initial_soc, current, voltage, Config, soc_values, ocv_values, time_data);
    
    % 비용 계산
    cost = sum(residuals.^2);
end

% Vt-Vest 잔차 계산 함수
function [residuals, Vt_est] = voltage_error(params, initial_soc, current, voltage, Config, soc_values, ocv_values, time_data)
    R0 = params(1);  % 고정된 R0
    R1 = params(2);  % 최적화된 R1
    C1 = params(3);  % 최적화된 C1

    num_samples = length(current);  % 샘플 개수
    SOC_est = zeros(num_samples, 1);  % SOC 추정값
    V1_est = zeros(num_samples, 1);   % V1 추정값
    Vt_est = zeros(num_samples, 1);   % 단자 전압 추정값
    residuals = zeros(num_samples, 1); % 잔차 벡터 초기화

    SOC_est(1) = initial_soc;  % 초기 SOC 설정
    dt_1 = time_data(2) - time_data(1);  % 첫 번째 샘플의 시간 간격
    V1_est(1) = current(1) * R1 * (1 - exp(-dt_1 / (R1 * C1)));  % 초기 V1 추정
    Vt_est(1) = voltage(1);  % 초기 전압 설정
    residuals(1) = Vt_est(1) - voltage(1);  % 초기 잔차 계산

    for k = 2:num_samples
        ik = current(k);
        dt_k = time_data(k) - time_data(k-1);  % 현재 샘플과 이전 샘플 간의 시간 차이

        SOC_est(k) = SOC_est(k-1) + (dt_k / (Config.cap * 3600)) * Config.coulomb_efficient * ik;
        V1_est(k) = exp(-dt_k / (R1 * C1)) * V1_est(k-1) + (1 - exp(-dt_k / (R1 * C1))) * ik * R1;
        Vt_est(k) = interp1(soc_values, ocv_values, SOC_est(k), 'linear', 'extrap') + V1_est(k) + R0 * ik;
        
        residuals(k) = Vt_est(k) - voltage(k);  % 각 샘플에 대한 잔차 계산
    end
end

