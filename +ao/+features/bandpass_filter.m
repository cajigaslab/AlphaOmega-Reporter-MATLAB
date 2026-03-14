function filtered = bandpass_filter(values, fs_hz, hp_hz, lp_hz, order)
%BANDPASS_FILTER  Zero-phase Butterworth bandpass filter.
%
%   filtered = ao.features.bandpass_filter(values, fs_hz, 300, 6000, 3)

    arguments
        values (:,1) double
        fs_hz  (1,1) double
        hp_hz  (1,1) double
        lp_hz  (1,1) double
        order  (1,1) double = 3
    end

    nyq = fs_hz / 2;
    [b, a] = butter(order, [hp_hz/nyq, lp_hz/nyq], 'bandpass');
    filtered = filtfilt(b, a, values);
end
