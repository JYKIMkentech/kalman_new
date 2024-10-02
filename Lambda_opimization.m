clc; clear; close all;

% Parameters
n = 21;  % Number of RC elements
t = 0:0.01:100;  % Time vector
dt = t(2) - t(1);

% Synthetic current data (sum of sine waves)
A1 = 1; 
A2 = 1;
A3 = 1;
T1 = 1;
T2 = 5;
T3 = 20;
I1 = A1 * sin(2 * pi * t / T1);  
I2 = A2 * sin(2 * pi * t / T2);
I3 = A3 * sin(2 * pi * t / T3);
ik = I1 + I2 + I3;  % Total current

% Parameters for the true DRT (R_discrete)
mu = 10;
sigma = 5;
tau_discrete = linspace(0.01, 20, n);  % Discrete tau values

% Calculate true R_discrete using a normal distribution
R_discrete = normpdf(tau_discrete, mu, sigma);
R_discrete = R_discrete / max(R_discrete);  % Normalize to max value of 1

% Initialize voltage
V_est = zeros(1, length(t));  % Estimated voltage
R0 = 0.1;  % Internal resistance
OCV = 0;   % Open circuit voltage
V_RC = zeros(n, length(t));  % RC voltages for each element

%% Initial voltage calculation (first time step)
for i = 1:n
    V_RC(i, 1) = ik(1) * R_discrete(i) * (1 - exp(-dt / tau_discrete(i))); % ik(1) = 0 
end
V_est(1) = OCV + R0 * ik(1) + sum(V_RC(:, 1));

%% Discrete-time voltage calculation for subsequent time steps
for k = 2:length(t)
    for i = 1:n
        % Calculate RC voltages based on previous time step
        V_RC(i, k) = exp(-dt / tau_discrete(i)) * V_RC(i, k-1) + R_discrete(i) * (1 - exp(-dt / tau_discrete(i))) * ik(k);       
    end
    V_est(k) = OCV + R0 * ik(k) + sum(V_RC(:, k));
end

%% Add noise to the voltage
rng(0);  % 시드 설정 (재현성을 위해)
noise_level = 0.01;
Noise = noise_level * randn(size(V_est));
V_noisy = V_est + Noise;

%% Regularization parameters and L matrix
% Construct the first derivative matrix L
L = zeros(n-1, n);
for i = 1:n-1
    L(i, i) = -1;
    L(i, i+1) = 1;
end

%% 교차 검증을 통한 람다 값 최적화
% 테스트할 람다 값 설정
lambda_values = [0, 0.001, 0.01, 0.1, 1];
num_lambdas = length(lambda_values);

% 시계열 교차 검증 설정
N = length(t);
initial_train_ratio = 0.5;  % 초기 훈련 데이터 비율
initial_train_size = floor(N * initial_train_ratio);
min_val_size = 50;  % 최소 검증 세트 크기
split_indices = initial_train_size:min_val_size:N - min_val_size;

% 람다 값별 평균 검증 오류 저장
mean_val_errors = zeros(num_lambdas, 1);

for l = 1:num_lambdas
    lambda = lambda_values(l);
    val_errors = zeros(length(split_indices), 1);
    
    for s = 1:length(split_indices)
        train_end = split_indices(s);
        val_start = train_end + 1;
        val_end = min(val_start + min_val_size - 1, N);
        
        % 훈련 및 검증 데이터 분할
        ik_train = ik(1:train_end);
        V_noisy_train = V_noisy(1:train_end);
        t_train = t(1:train_end);
        
        ik_val = ik(val_start:val_end);
        V_noisy_val = V_noisy(val_start:val_end);
        t_val = t(val_start:val_end);
        
        % W 행렬 생성 (훈련 데이터용)
        W_train = construct_W_matrix(tau_discrete, ik_train, dt);
        
        % y 벡터 생성 (훈련 데이터용)
        y_train = V_noisy_train' - OCV - R0 * ik_train';
        
        % 정규 방정식을 사용하여 R 추정
        A = W_train' * W_train + lambda * L' * L;
        b = W_train' * y_train;
        R_est = A \ b;
        
        % 음수 값 제거
        R_est(R_est < 0) = 0;
        
        % 검증 세트에 대한 예측 전압 계산
        V_est_val = calculate_voltage(R_est', tau_discrete, ik_val, dt, n, R0, OCV, t_val);
        
        % 검증 오류 계산 (MSE)
        val_errors(s) = mean((V_noisy_val - V_est_val).^2);
    end
    
    % 평균 검증 오류 저장
    mean_val_errors(l) = mean(val_errors);
    fprintf('Lambda = %f, Mean Validation Error = %f\n', lambda, mean_val_errors(l));
end

% 최적의 람다 값 선택
[~, optimal_idx] = min(mean_val_errors);
optimal_lambda = lambda_values(optimal_idx);
fprintf('Optimal Lambda: %f\n', optimal_lambda);

%% 최적의 람다로 전체 데이터에 대해 모델 학습
lambda = optimal_lambda;

% W 행렬 생성 (전체 데이터용)
W = construct_W_matrix(tau_discrete, ik, dt);

% y 벡터 생성 (전체 데이터용)
y = V_noisy' - OCV - R0 * ik';

% 정규 방정식을 사용하여 R 추정
A = W' * W + lambda * L' * L;
b = W' * y;
R_optimal = A \ b;
R_optimal(R_optimal < 0) = 0;

%% 결과 비교 및 시각화
figure;
hold on;

plot(tau_discrete, R_discrete, 'b-', 'LineWidth', 2);  % 실제 DRT (파란색)
stem(tau_discrete, R_discrete, 'bo', 'LineWidth', 1.5);

plot(tau_discrete, R_optimal, 'r-', 'LineWidth', 2);  % 최적 람다를 사용한 DRT (빨간색)
stem(tau_discrete, R_optimal, 'ro', 'LineWidth', 1.5);

xlabel('\tau (Time Constant)');
ylabel('R (Resistance)');
legend('True DRT', 'True DRT Points', 'Optimized DRT', 'Optimized DRT Points');
title(['Optimal Lambda: ', num2str(optimal_lambda)]);
grid on;
hold off;

%% Functions

% W 행렬을 생성하는 함수
function W = construct_W_matrix(tau_discrete, ik, dt)
    n = length(tau_discrete);
    N = length(ik);
    W = zeros(N, n);
    for i = 1:n
        for k = 1:N
            if k == 1
                W(k, i) = ik(k) * (1 - exp(-dt / tau_discrete(i)));
            else
                W(k, i) = exp(-dt / tau_discrete(i)) * W(k-1, i) + ik(k) * (1 - exp(-dt / tau_discrete(i)));
            end
        end
    end
end

% 주어진 R_discrete로 V_est를 계산하는 함수
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
