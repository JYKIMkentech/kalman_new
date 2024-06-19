function v_ocv = ocv_soc(soc)
    % OCV-SOC relationship
    v_ocv = 2.58 * soc + 3.81 * exp(-35.8 * soc) - 0.3 * exp(-9.8 * soc);
end
