function filtered = apply_notch_filter(values, fs_hz, center_hz, width_hz)
%APPLY_NOTCH_FILTER  Zero-phase IIR notch filter.
%
%   filtered = ao.features.apply_notch_filter(values, fs_hz, 60, 2)

    arguments
        values    (:,1) double
        fs_hz     (1,1) double
        center_hz (1,1) double
        width_hz  (1,1) double = 2.0
    end

    Q = center_hz / width_hz;
    wo = center_hz / (fs_hz / 2);  % normalized frequency

    [b, a] = iirnotch(wo, wo/Q);
    filtered = filtfilt(b, a, values);
end
