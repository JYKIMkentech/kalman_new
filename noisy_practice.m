% 파라미터 설정
n = 21;                    % 상태의 개수
sigma = 1e-3;              % 1 시그마에 5개 구간을 포함하도록 설정 (1 mA)
noise_levels = linspace(-2e-3, 2e-3, n); % -2mA에서 +2mA로 20개 구간

% 전이 상태 행렬 P 초기화
P = zeros(n, n);

% 전이 확률 계산
for i = 1:n
    for j = 1:n
        P(i, j) = normpdf(noise_levels(j), noise_levels(i), sigma);
    end
    % 각 행의 합이 1이 되도록 정규화
    P(i, :) = P(i, :) / sum(P(i, :));
end

% 노이즈 벡터 Q 초기화
Q = noise_levels';

% True current 데이터 설정 (사용자가 제공한 값)
true_current = [0 0 0 5 5 5 0 0 0 -3 -3 -3] * 1e-3; % A 단위로 스케일링 (mA -> A)

% 배터리 용량 및 초기 SOC 설정
battery_capacity = 2.9; % 배터리 용량 2.9 Ah
initial_soc = 1;        % 초기 SOC는 100% (1.0)

% 시뮬레이션 스텝 수 및 시간 간격 설정
num_steps = length(true_current);
dt = 1; % 시간 간격 (단위: 시간, 예를 들어 1초)

% SOC 계산을 위한 적분 변수 초기화
true_soc = zeros(1, num_steps);
noisy_current = zeros(1, num_steps);
noisy_soc = zeros(1, num_steps);

% 초기 SOC 설정
true_soc(1) = initial_soc;
noisy_soc(1) = initial_soc;

% 마르코브 체인을 이용한 noisy current 생성 및 SOC 계산
for t = 2:num_steps
    % 상태 전이
    current_state = randsample(1:n, 1, true, P(current_state, :));
    
    % 노이즈 추가 (mA를 A로 스케일링하여 추가)
    noisy_current(t) = true_current(t) + Q(current_state);
    
    % True SOC 계산 (전류를 적분하여 계산)
    true_soc(t) = true_soc(t-1) - (true_current(t) * dt) / battery_capacity;
    
    % Noisy SOC 계산 (전류를 적분하여 계산)
    noisy_soc(t) = noisy_soc(t-1) - (noisy_current(t) * dt) / battery_capacity;
end

% 결과 시각화
figure;
subplot(2,1,1);
plot(1:num_steps, true_soc * 100, 'b-', 'LineWidth', 2); % True SOC (%)
hold on;
plot(1:num_steps, noisy_soc * 100, 'r--', 'LineWidth', 2); % Noisy SOC (%)
xlabel('Time Step');
ylabel('SOC (%)');
title('True SOC vs Noisy SOC');
legend('True SOC', 'Noisy SOC');
grid on;
hold off;

subplot(2,1,2);
plot(1:num_steps, true_current * 1e3, 'b-', 'LineWidth', 2); % True current (mA)
hold on;
plot(1:num_steps, noisy_current * 1e3, 'r--', 'LineWidth', 2); % Noisy current (mA)
xlabel('Time Step');
ylabel('Current (mA)');
title('True Current vs Noisy Current');
legend('True Current', 'Noisy Current');
grid on;
hold off;

