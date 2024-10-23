clear; clc; close all;

%% 1. UDDS 주행 데이터 로드
% UDDS 주행 데이터를 로드합니다.
load('C:\Users\deu04\OneDrive\바탕 화면\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat');
udds_current = meas.Current;  % 전류 데이터 (A)
udds_voltage = meas.Voltage;  % 전압 데이터 (V)
udds_time = meas.Time;        % 시간 데이터 (s)

% 시간 벡터가 0에서 시작하고 연속적이도록 보정합니다.
udds_time = udds_time - udds_time(1);

%% 2. SOC-OCV 데이터 로드
% SOC-OCV 데이터를 로드합니다.
load('soc_ocv.mat', 'soc_ocv');
soc_values = soc_ocv(:, 1);  % SOC 값
ocv_values = soc_ocv(:, 2);  % OCV 값

%% 3. SOC 계산
% 초기 SOC 설정
SOC_initial = 0.9901;  % 초기 SOC를 설정합니다.

% 전류 데이터의 시간 간격(delta_t) 계산
delta_t = [0; diff(udds_time)];  % 각 측정 사이의 시간 간격

% 배터리 용량 (Ah를 Coulomb로 변환)
Q_battery = 2.9 * 3600;  % 배터리 용량 (2.9Ah)

% 시간에 따른 SOC 계산
udds_SOC = SOC_initial + cumsum(udds_current .* delta_t) / Q_battery;

%% 4. UDDS 주기 시작점과 끝점 탐지 (Sequential Approach)
% 초기 설정
initial_cycle_start_time = 0;          % 첫 주기 시작 시간 (0초)
next_cycle_duration = 1370;            % 다음 주기까지의 예상 시간 (초)
current_at_zero = udds_current(1);     % 시간 0에서의 전류
current_threshold = 0.001;             % 전류 유사성 임계값 (A)
time_window = 10;                      % 목표 시간 주변의 시간 창 (초)
total_time = udds_time(end);           % 총 시간

% 주기 시작 인덱스를 저장할 배열 초기화 (첫 주기 시작점 포함)
cycle_start_indices = 1; % 첫 주기 시작점 (인덱스 1)

% 현재 탐색할 주기 시작 시간
current_cycle_start_time = initial_cycle_start_time;

% 주기 탐색 루프
while true
    % 다음 주기 시작 시간 예상
    t_target = current_cycle_start_time + next_cycle_duration;

    % t_target이 총 시간을 초과하면 루프 종료
    if t_target > total_time
        break;
    end

    % 시간 창 내의 인덱스 찾기
    indices_in_window = find(abs(udds_time - t_target) <= time_window);

    % 해당 인덱스 중 전류가 시간 0의 전류와 유사한 인덱스 찾기
    indices_with_similar_current = indices_in_window(abs(udds_current(indices_in_window) - current_at_zero) <= current_threshold);

    if ~isempty(indices_with_similar_current)
        % 목표 시간에 가장 가까운 인덱스 선택
        [~, idx_closest_time] = min(abs(udds_time(indices_with_similar_current) - t_target));
        index = indices_with_similar_current(idx_closest_time);

        % 인덱스 저장
        cycle_start_indices = [cycle_start_indices; index];

        % 다음 탐색을 위해 현재 주기 시작 시간 업데이트
        current_cycle_start_time = udds_time(index);
    else
        % 유사한 전류를 찾지 못한 경우, 다음 주기 탐색을 종료
        break;
    end
end

% 주기 시작점을 중복 없이 정렬
cycle_start_indices = unique(cycle_start_indices, 'sorted');

%% 5. udds_data 구조체 배열 생성
% 각 주기의 데이터를 동일한 필드 이름을 가진 구조체 배열로 저장
num_cycles = length(cycle_start_indices);

% udds_data 구조체 배열 초기화
udds_data = struct('V', {}, 'I', {}, 't', {}, 'SOC', {});

for i = 1:num_cycles
    if i < num_cycles
        start_idx = cycle_start_indices(i);
        end_idx = cycle_start_indices(i+1) - 1;
    else
        start_idx = cycle_start_indices(i);
        end_idx = length(udds_time);
    end

    % 각 주기의 데이터를 구조체 배열의 요소로 저장
    udds_data(i).V = udds_voltage(start_idx:end_idx);
    udds_data(i).I = udds_current(start_idx:end_idx);
    udds_data(i).Time_duration = udds_time(start_idx:end_idx);
    udds_data(i).t = udds_time(start_idx:end_idx) - udds_time(start_idx);
    udds_data(i).SOC = udds_SOC(start_idx:end_idx);
end

%% 6. 찾은 주기 시작 시간을 fprintf로 출력
fprintf('Detected Cycle Start Times:\n');
fprintf('--------------------------\n');
for i = 1:length(cycle_start_indices)
    cycle_num = i; % 주기 번호 (Cycle 1, Cycle 2, ...)
    time_sec = udds_time(cycle_start_indices(i)); % 주기 시작 시간 (초)

    % 마지막 주기의 경우 'Cycle N.xx'와 같이 특별한 이름을 지정할 수 있습니다.
    if i == length(cycle_start_indices)
        fprintf('Cycle %d start time: %.2f s (Cycle %.2f)\n', cycle_num, time_sec, cycle_num - 1 + 0.38);
    else
        fprintf('Cycle %d start time: %.2f s\n', cycle_num, time_sec);
    end
end

%% 7. udds_data 구조체 내용 확인
% 구조체 배열의 각 요소에 접근하여 데이터를 확인할 수 있습니다.
disp('udds_data 구조체 배열의 내용:');
disp('---------------------------');

for i = 1:num_cycles
    fprintf('Cycle %d:\n', i);
    fprintf('  Time length: %d samples\n', length(udds_data(i).t));
    fprintf('  Voltage length: %d samples\n', length(udds_data(i).V));
    fprintf('  Current length: %d samples\n', length(udds_data(i).I));
    fprintf('  SOC length: %d samples\n', length(udds_data(i).SOC));
end

%% 8. Plot with Correct Legend

figure;
hold on;

% Plot Current on the left y-axis
yyaxis left
h_current = plot(udds_time, udds_current, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Current (A)');
ylabel('Current (A)');
ylim([min(udds_current)-0.5, max(udds_current)+0.5]);

% Plot SOC on the right y-axis
yyaxis right
h_soc = plot(udds_time, udds_SOC, 'g-', 'LineWidth', 1.5, 'DisplayName', 'SOC');
ylabel('SOC');
ylim([min(udds_SOC)-0.05, max(udds_SOC)+0.05]);

% Initialize arrays to store handles for cycle markers and lines
h_cycle_markers = [];
h_cycle_lines = [];

% Plot cycle boundaries
for i = 1:length(cycle_start_indices)
    x = udds_time(cycle_start_indices(i));
    y_current = udds_current(cycle_start_indices(i));
    y_soc = udds_SOC(cycle_start_indices(i));

    % Plot 'ro' marker on left y-axis for Current
    yyaxis left
    h_marker_current = plot(x, y_current, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    h_cycle_markers(end+1) = h_marker_current; 

    % Plot 'ro' marker on right y-axis for SOC
    yyaxis right
    h_marker_soc = plot(x, y_soc, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    h_cycle_markers(end+1) = h_marker_soc; 

    % Plot vertical dashed line on left y-axis
    yyaxis left
    h_line = plot([x x], ylim, 'k--', 'LineWidth', 1);
    h_cycle_lines(end+1) = h_line; 

    % Add cycle number labels
    if i < length(cycle_start_indices)
        midpoint = (udds_time(cycle_start_indices(i)) + udds_time(cycle_start_indices(i+1))) / 2;
    else
        midpoint = udds_time(cycle_start_indices(i)) + (udds_time(end) - udds_time(cycle_start_indices(i))) / 2;
    end

    % Adjust label for the last cycle
    if i == length(cycle_start_indices)
        label = sprintf('Cycle %.2f', i - 1 + 0.38);
    else
        label = sprintf('Cycle %d', i);
    end

    % Add text label above the plot
    text(midpoint, max(udds_current)+0.3, label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

% Reset to left y-axis for consistency
yyaxis left

% Create dummy plots for legend entries not directly tied to a single plot
h_cycle_marker_dummy = plot(NaN, NaN, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Cycle Boundaries');
h_cycle_line_dummy = plot(NaN, NaN, 'k--', 'LineWidth', 1, 'DisplayName', 'Cycle Lines');

% Combine all legend handles
legend_handles = [h_current, h_soc, h_cycle_marker_dummy, h_cycle_line_dummy];

% Add legend to the plot
legend(legend_handles, 'Location', 'best');

% Finalize plot settings
xlabel('Time (s)');
title('UDDS Current and SOC Profile with Cycle Boundaries');
grid on;
hold off;


%% Save

save('udds_data.mat', 'udds_data');

