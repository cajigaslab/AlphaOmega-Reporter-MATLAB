function feat = compute_features(values, fs_hz)
%COMPUTE_FEATURES  Compute all time-domain and spectral features for a signal.
%
%   feat = ao.features.compute_features(values, fs_hz)
%
%   Returns a struct with fields:
%       .mer_rms            - root mean square (µV)
%       .kurtosis           - excess kurtosis (0 for Gaussian)
%       .spectral_entropy   - normalized spectral entropy [0, 1]
%       .zero_crossing_rate - zero crossings per second (Hz)
%       .line_length        - mean absolute first difference
%       .hjorth_activity    - signal variance
%       .hjorth_mobility    - std(diff)/std(signal)
%       .hjorth_complexity  - mobility(diff)/mobility(signal)

    arguments
        values (:,1) double
        fs_hz  (1,1) double
    end

    v = values(isfinite(values));
    N = numel(v);

    % --- RMS ---
    feat.mer_rms = sqrt(mean(v.^2));

    % --- Excess kurtosis ---
    if N > 4
        mu = mean(v);
        sigma = std(v);
        if sigma > 0
            feat.kurtosis = mean((v - mu).^4) / sigma^4 - 3;
        else
            feat.kurtosis = 0;
        end
    else
        feat.kurtosis = 0;
    end

    % --- Spectral entropy ---
    if N >= 128
        nperseg = min(512, 2^floor(log2(N)));
        [pxx, ~] = pwelch(v, nperseg, floor(nperseg/2), [], fs_hz, 'onesided');
        p = pxx / sum(pxx);
        p = p(p > 0);
        feat.spectral_entropy = -sum(p .* log(p)) / log(numel(p));
    else
        feat.spectral_entropy = 1.0;
    end

    % --- Zero crossing rate ---
    if N > 1
        centered = v - mean(v);
        crossings = sum(abs(diff(sign(centered))) > 0);
        duration_s = N / fs_hz;
        feat.zero_crossing_rate = crossings / duration_s;
    else
        feat.zero_crossing_rate = 0;
    end

    % --- Line length ---
    if N > 1
        feat.line_length = mean(abs(diff(v)));
    else
        feat.line_length = 0;
    end

    % --- Hjorth parameters ---
    if N > 2
        d1 = diff(v);
        d2 = diff(d1);
        feat.hjorth_activity   = var(v);
        mob_signal = std(d1) / std(v);
        feat.hjorth_mobility   = mob_signal;
        if mob_signal > 0
            mob_d1 = std(d2) / std(d1);
            feat.hjorth_complexity = mob_d1 / mob_signal;
        else
            feat.hjorth_complexity = 0;
        end
    else
        feat.hjorth_activity   = 0;
        feat.hjorth_mobility   = 0;
        feat.hjorth_complexity = 0;
    end
end
