function [SOC_est, V1_est, Vt_est, P] = kalman_filter(SOC_est, V1_est, Vt_true, ik, Config, P)
    Q = [1e-5 0; 0 1e-5]; % Process noise covariance
    R = 1e-2; % Measurement noise covariance
    
    % Prediction step
    SOC_pred = SOC_est - (Config.dt / Config.cap) * Config.coulomb_efficient * ik;
    V1_pred = exp(-Config.dt / (Config.R1 * Config.C1)) * V1_est + (1 - exp(-Config.dt / (Config.R1 * Config.C1))) * ik * Config.R1;
    X_pred = [SOC_pred; V1_pred];
    
    % State transition matrix
    A = [1, 0; 0, exp(-Config.dt / (Config.R1 * Config.C1))];
    
    % Process covariance update
    P = A * P * A' + Q;
    
    % Measurement prediction
    Vt_pred = OCV(SOC_pred) - Config.R1 * V1_pred - Config.R0 * ik;
    
    % Measurement update
    K = P * [1; 0] / (P(1, 1) + R); % Kalman gain
    z = Vt_true - Vt_pred; % Measurement residual
    X_pred = X_pred + K * z;
    P = (eye(2) - K * [1, 0]) * P;
    
    % Save estimates
    SOC_est = X_pred(1);
    V1_est = X_pred(2);
    Vt_est = Vt_pred + K(1) * z; % Estimated Vt
end

function v_ocv = OCV(soc)
    % OCV-SOC relationship
    v_ocv = 2.58 * soc + 3.81 * exp(-35.8 * soc) - 0.3 * exp(-9.8 * soc);
end
