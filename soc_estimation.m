function [SOC_est, V1_est, Vt_est] = soc_estimation(SOC_est, V1_est, Vt_true, ik, Config, P)
    [SOC_est, V1_est, Vt_est, P] = kalman_filter(SOC_est, V1_est, Vt_true, ik, Config, P);
end
