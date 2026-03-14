function result = compute_psd(values, fs_hz, cfg)
%COMPUTE_PSD  Power spectral density via multitaper or Welch method.
%
%   result = ao.features.compute_psd(values, fs_hz, cfg)
%
%   Returns:
%       result.freq_hz  - frequency vector (column)
%       result.psd      - power spectral density (column)
%       result.psd_db   - 10*log10(psd) in dB (column)

    arguments
        values  (:,1) double
        fs_hz   (1,1) double
        cfg     struct
    end

    % Strip non-finite
    values = values(isfinite(values));
    if numel(values) < 64
        result.freq_hz = [];
        result.psd     = [];
        result.psd_db  = [];
        return;
    end

    fmin = cfg.lfp.fmin_hz;
    fmax = cfg.lfp.fmax_hz;

    switch lower(cfg.lfp.spectrum_method)
        case 'multitaper'
            result = multitaper_psd(values, fs_hz, fmin, fmax, cfg.lfp.multitaper_nw);
        case 'welch'
            result = welch_psd(values, fs_hz, fmin, fmax, ...
                cfg.lfp.welch_nperseg, cfg.lfp.welch_noverlap);
        otherwise
            error('ao:features:compute_psd', ...
                'Unknown spectrum method: %s', cfg.lfp.spectrum_method);
    end
end

%% ========================================================================
function result = multitaper_psd(values, fs_hz, fmin, fmax, nw)
    N = numel(values);
    nfft = max(2^nextpow2(N), 256);
    kmax = max(floor(2*nw - 1), 1);

    % DPSS tapers (requires Signal Processing Toolbox)
    [tapers, ~] = dpss(N, nw, kmax);

    % Compute PSD for each taper and average
    psd_sum = zeros(nfft, 1);
    for k = 1:kmax
        tapered = values .* tapers(:, k);
        X = fft(tapered, nfft);
        taper_energy = sum(tapers(:, k).^2);
        power = (abs(X).^2) / (fs_hz * taper_energy);
        psd_sum = psd_sum + power;
    end
    psd_full = psd_sum / kmax;

    % One-sided spectrum
    n_onesided = floor(nfft/2) + 1;
    freq = (0:n_onesided-1)' * (fs_hz / nfft);
    psd = psd_full(1:n_onesided);
    psd(2:end-1) = 2 * psd(2:end-1);  % double non-DC, non-Nyquist bins

    % Frequency mask
    mask = freq >= fmin & freq <= fmax;
    result.freq_hz = freq(mask);
    result.psd     = psd(mask);
    result.psd_db  = 10 * log10(psd(mask) + eps);
end

%% ========================================================================
function result = welch_psd(values, fs_hz, fmin, fmax, nperseg, noverlap)
    N = numel(values);
    nperseg = min(nperseg, 2^floor(log2(N)));
    nperseg = max(nperseg, 128);
    noverlap = min(noverlap, nperseg - 1);

    [psd, freq] = pwelch(values, nperseg, noverlap, [], fs_hz, 'onesided');

    mask = freq >= fmin & freq <= fmax;
    result.freq_hz = freq(mask);
    result.psd     = psd(mask);
    result.psd_db  = 10 * log10(psd(mask) + eps);
end
