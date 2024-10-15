clc; clear; close all;

% tau는 로그 정규 분포를 따름
mean_lognormal = 10;
std_lognormal = 2;

% ln(tau)는 정규분포를 따름
sigma_squared = log(1 + (std_lognormal^2) / (mean_lognormal^2));
sigma = sqrt(sigma_squared); % ln(tau)의 표준편차
mu = log(mean_lognormal) - (sigma_squared / 2); % ln(tau)의 평균

% theta = ln(tau) 
ln_tau_min = log(1);  % theta_min = log(0.1) 
ln_tau_max = log(20);   % theta_max = log(20) 
ln_tau_values = linspace(ln_tau_min, ln_tau_max, 100); % ln(tau)을 일정 간격으로 나눔

% ln(tau) = theta 를 tau로 변환하여 tau = exp(theta) 
tau_values = exp(ln_tau_values); % exp(theta)

% 로그정규 분포 PDF 계산  (tau 기준)
lognormal_pdf = lognpdf(tau_values, mu, sigma);

% 정규분포 PDF 계산 (ln(tau) 기준 = theta 기준)
normal_pdf = normpdf(ln_tau_values, mu, sigma);

% Plot
figure;

% tau의 log normal 분포 
subplot(2, 1, 1);
plot(tau_values, lognormal_pdf, 'r-', 'LineWidth', 2);
xlabel('\tau');
ylabel('resistance');
title('Lognormal Distribution of \tau (DRT)');
set(gca, 'XScale', 'log'); % tau 축을 로그 스케일로 설정
grid on;

% ln(tau)의 normal 분포 
subplot(2, 1, 2);
plot(ln_tau_values, normal_pdf, 'b-', 'LineWidth', 2);
xlabel('ln(\tau)');
ylabel('Gamma');
title('Normal Distribution of ln(\tau) (Reference DRT)');
grid on;

