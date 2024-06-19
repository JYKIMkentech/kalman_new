function [SOC_true, Vt_true] = battery_model(SOC_prev, V1_prev, ik, Config)
    SOC_true = SOC_prev - (Config.dt / Config.cap) * Config.coulomb_efficient * ik;
    Vt_true = OCV(SOC_true) - Config.R1 * V1_prev - Config.R0 * ik;
end

function v_ocv = OCV(soc)
    % OCV-SOC relationship
    v_ocv = 2.58 * soc + 3.81 * exp(-35.8 * soc) - 0.3 * exp(-9.8 * soc);
end
