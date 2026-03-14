function result = detect_spikes(values, fs_hz, cfg)
%DETECT_SPIKES  Threshold-based spike detection on MER signal.
%
%   result = ao.detection.detect_spikes(values, fs_hz, cfg)
%
%   Algorithm:
%     1. Bandpass filter (default 300–6000 Hz, 3rd-order Butterworth)
%     2. Estimate noise: sigma = MAD / 0.6745
%     3. Threshold = threshold_mad × sigma
%     4. Find threshold crossings, peak within ±0.5 ms
%     5. Enforce dead time between events
%
%   Returns:
%       result.spike_times_s    - spike times in seconds (column)
%       result.filtered         - bandpass-filtered signal (column)
%       result.sigma_uv         - noise estimate (µV)
%       result.threshold_uv     - threshold used (µV)
%       result.spikeband_rms    - RMS of filtered signal
%       result.n_spikes         - number of detected spikes

    arguments
        values (:,1) double
        fs_hz  (1,1) double
        cfg    struct
    end

    hp_hz     = cfg.detect.filter_hp_hz;
    lp_hz     = cfg.detect.filter_lp_hz;
    order     = cfg.detect.filter_order;
    thr_mad   = cfg.detect.threshold_mad;
    peak_sign = cfg.detect.peak_sign;
    dead_ms   = cfg.detect.dead_time_ms;

    % --- Bandpass filter ---
    filtered = ao.features.bandpass_filter(values, fs_hz, hp_hz, lp_hz, order);

    % --- Noise estimate (MAD) ---
    sigma = median(abs(filtered - median(filtered))) / 0.6745;
    threshold_uv = thr_mad * sigma;

    % --- Find threshold crossings ---
    search_radius = round(0.5e-3 * fs_hz);  % ±0.5 ms
    dead_samples  = round(dead_ms / 1000 * fs_hz);

    switch lower(peak_sign)
        case 'neg'
            exceeds = filtered < -threshold_uv;
        case 'pos'
            exceeds = filtered > threshold_uv;
        case 'both'
            exceeds = abs(filtered) > threshold_uv;
        otherwise
            exceeds = filtered < -threshold_uv;
    end

    % Group contiguous crossings into events, find peak within each
    spike_indices = find_peaks_in_crossings(filtered, exceeds, ...
        search_radius, peak_sign);

    % --- Enforce dead time ---
    if numel(spike_indices) > 1
        spike_indices = enforce_dead_time(spike_indices, filtered, dead_samples);
    end

    % --- Build result ---
    result.spike_times_s  = (spike_indices - 1) / fs_hz;  % 0-indexed time
    result.filtered       = filtered;
    result.sigma_uv       = sigma;
    result.threshold_uv   = threshold_uv;
    result.spikeband_rms  = sqrt(mean(filtered.^2));
    result.n_spikes       = numel(spike_indices);
end

%% ========================================================================
function peaks = find_peaks_in_crossings(filtered, exceeds, radius, peak_sign)
    N = numel(filtered);
    % Label contiguous crossing groups
    d = diff([0; exceeds(:); 0]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    peaks = zeros(numel(starts), 1);
    for k = 1:numel(starts)
        % Extend search window by radius
        lo = max(1, starts(k) - radius);
        hi = min(N, stops(k) + radius);
        segment = filtered(lo:hi);

        switch lower(peak_sign)
            case 'neg'
                [~, idx] = min(segment);
            case 'pos'
                [~, idx] = max(segment);
            case 'both'
                [~, idx] = max(abs(segment));
            otherwise
                [~, idx] = min(segment);
        end
        peaks(k) = lo + idx - 1;
    end

    % Remove duplicates
    peaks = unique(peaks);
end

%% ========================================================================
function kept = enforce_dead_time(indices, filtered, dead_samples)
    indices = sort(indices);
    kept = indices(1);
    for k = 2:numel(indices)
        if indices(k) - kept(end) >= dead_samples
            kept(end+1) = indices(k); %#ok<AGROW>
        else
            % Keep the one with larger amplitude
            if abs(filtered(indices(k))) > abs(filtered(kept(end)))
                kept(end) = indices(k);
            end
        end
    end
end
