clc;clear;close all;

%% 1. Define Time Constants on Logarithmic Scale (Natural Log)
a = 0.01;        % 최소 tau
b = 100;         % 최대 tau
n = 50;          % 타우 개수 (필요에 따라 변경)
theta_j = linspace(log(a), log(b), n); % 로그 도메인에서 균등 간격
tau_discrete_log = exp(theta_j);        % tau_j 계산

%% 2. True DRT Parameters (Gamma) - Lognormal 분포 사용
% 로그정규분포의 평균과 표준편차 설정
mean_tau = 10;      % tau의 평균
std_tau = 5;        % tau의 표준편차

% 기저 정규분포의 파라미터 mu와 sigma 계산
sigma_log = sqrt(log(1 + (std_tau / mean_tau)^2));
mu_log = log(mean_tau) - (sigma_log^2) / 2;

%% 3. Lognormal PDF 계산
% tau_discrete_log에 대한 Gamma_true 계산 (로그정규분포의 PDF)
Gamma_true = lognpdf(tau_discrete_log, mu_log, sigma_log);

%% 4. Plotting the result (x축은 ln(tau))
figure;
plot(log(tau_discrete_log), Gamma_true, 'b-', 'LineWidth', 2);
xlabel('ln(\tau)');
ylabel('Gamma(\tau)');
title('Lognormal Distribution of Time Constants (DRT)');
grid on;
