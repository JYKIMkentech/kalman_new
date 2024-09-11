% 데이터 초기화
t = [0, 1, 2, 3, 4, 5]; % 시간 (초)
V_t = [2.0, 2.1, 2.2, 2.1, 2.0, 1.9]; % 전압 (볼트)
I_t = [0.5, 0.52, 0.54, 0.52, 0.5, 0.48]; % 전류 (암페어)

% 임피던스 계산 (Z = V / I)
Z_t = V_t ./ I_t;

% 임피던스 값 출력
disp('시간 도메인 임피던스:');
disp(Z_t);

% FFT를 사용하여 주파수 도메인으로 변환
Z_f = fft(Z_t);

% FFT 결과 출력
disp('FFT 결과 (복소수 임피던스):');
disp(Z_f');

% 주파수 벡터 계산 (기본적인 설정으로 간단히 수행)
N = length(Z_t);  % 데이터 길이
Fs = 1 / (t(2) - t(1));  % 샘플링 주파수
frequencies = (0:N-1) * (Fs / N);  % 주파수 벡터

% 주파수 도메인에서의 임피던스 출력
disp('주파수 (Hz):');
disp(frequencies);
disp('주파수 도메인 임피던스:');
disp(Z_f);

% DRT 분석을 위한 시간 상수 범위 설정
tau_min = 1e-3;  % 최소 시간 상수
tau_max = 1e1;   % 최대 시간 상수
n_tau = 50;     % 시간 상수 개수
tau = logspace(log10(tau_min), log10(tau_max), n_tau);  % 로그 스케일의 시간 상수

% Tikhonov 정규화를 적용한 DRT 계산
lambda = 0.01;  % 정규화 파라미터
H = zeros(size(tau));  % DRT 함수 초기화
for i = 1:n_tau
    G = exp(-t / tau(i));
    H(i) = sum(real(Z_f) .* G) / (sum(G.^2) + lambda);
end

% 결과 플로팅
figure;
semilogx(tau, H, 'LineWidth', 2);
xlabel('Time Constant (\tau) [s]');
ylabel('DRT Distribution h(\tau)');
title('Distribution of Relaxation Times (DRT)');
grid on;
