clc;clear; close all;

% tau는 log normal 분포를 따름
mean_lognormal = 10;
std_lognormal = 5;

% ln(tau)는 정규분포를 따름
sigma_squared = log(1 + (std_lognormal^2) / (mean_lognormal^2));
sigma = sqrt(sigma_squared); % ln(tau)의 표준편차
mu = log(mean_lognormal) - (sigma_squared / 2); % ln(tau)의 평균

% tau 범위 설정
%tau_values = logspace(log10(mean_lognormal/10), log10(mean_lognormal*10), 100);
tau_values = linspace(0.1, 50, 100);

% 로그정규 분포 PDF 계산 
lognormal_pdf = lognpdf(tau_values, mu, sigma);

% ln(tau) 범위 설정
ln_tau_values = log(tau_values);

% 정규분포 PDF 계산
normal_pdf = normpdf(ln_tau_values, mu, sigma);

% Plot
figure;

% tau의 log normal 분포 
subplot(2, 1, 1);
plot(tau_values, lognormal_pdf, 'r-', 'LineWidth', 2);
xlabel('\tau');
ylabel('Resistance (Gamma)');
title('Lognormal Distribution of \tau (DRT)');
grid on;

% ln (tau)의 normal 분포 
subplot(2, 1, 2);
plot(ln_tau_values, normal_pdf, 'b-', 'LineWidth', 2);
xlabel('ln(\tau)');
ylabel('Gamma');
title('Normal Distribution of ln(\tau) (Reference DRT)');
grid on;

