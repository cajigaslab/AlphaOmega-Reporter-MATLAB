function is_imp = detect_impedance(rms, kurtosis_val, spectral_entropy, median_rms, cfg)
%DETECT_IMPEDANCE  Composite 2-of-3 impedance check detection.
%
%   is_imp = ao.qc.detect_impedance(rms, kurtosis, spectral_entropy, median_rms, cfg)
%
%   Criteria (any 2 of 3 triggers a flag):
%     1. RMS ratio > 10×  (extreme signal power)
%     2. Kurtosis  > 20   (artifact-like distribution)
%     3. Spectral entropy < 0.4  (narrow-band, non-neural)

    arguments
        rms              (1,1) double
        kurtosis_val     (1,1) double
        spectral_entropy (1,1) double
        median_rms       (1,1) double
        cfg              struct
    end

    criteria = 0;

    % RMS ratio
    if median_rms > 0 && (rms / median_rms) > cfg.impedance.rms_ratio_threshold
        criteria = criteria + 1;
    end

    % Kurtosis
    if kurtosis_val > cfg.impedance.kurtosis_threshold
        criteria = criteria + 1;
    end

    % Spectral entropy
    if spectral_entropy < cfg.impedance.spectral_entropy_threshold
        criteria = criteria + 1;
    end

    is_imp = criteria >= cfg.impedance.min_criteria;
end
