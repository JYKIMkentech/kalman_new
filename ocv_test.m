clc;clear;clc

% Load C/20 charge/discharge data
load('C:\Users\USER\Desktop\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat');
current = meas.Current; % 전류 데이터
voltage = meas.Voltage; % 전압 데이터
time = meas.Time; % 시간 데이터

% Extract discharge portion
discharge_idx = current < 0;
discharge_current = current(discharge_idx);
discharge_voltage = voltage(discharge_idx);
discharge_time = time(discharge_idx);

% Initialize SOC calculation
%initial_capacity = 4.32019; % Ah, 초기 용량

% Calculate total charge (total Q) using trapezoidal integration
total_q = trapz(discharge_time, discharge_current) / 3600; % Ah

% Calculate SOC using real-time integration
SOC = zeros(length(discharge_current), 1);
SOC(1) = 1; % 초기 SOC는 100%로 가정

for k = 2:length(discharge_current)
    delta_t = discharge_time(k) - discharge_time(k-1);
    SOC(k) = SOC(k-1) - (discharge_current(k) * delta_t) / (total_q * 3600);
end

% Ensure SOC is within bounds [0, 1]
SOC = max(0, min(1, SOC));

% Match OCV values
OCV = discharge_voltage;

% Create SOC-OCV structure
soc_ocv = [SOC OCV];

% Save SOC-OCV structure for later use
save('soc_ocv.mat', 'soc_ocv');