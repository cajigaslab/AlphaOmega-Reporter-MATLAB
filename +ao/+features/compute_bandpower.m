function [power_db, std_db] = compute_bandpower(freq_hz, psd, psd_db, low_hz, high_hz)
%COMPUTE_BANDPOWER  Integrated band power in dB and within-band std.
%
%   [power_db, std_db] = ao.features.compute_bandpower(freq, psd, psd_db, 13, 30)
%
%   power_db  = 10*log10(trapz(freq, psd) over [low_hz, high_hz])
%   std_db    = std(psd_db) within the band

    arguments
        freq_hz (:,1) double
        psd     (:,1) double
        psd_db  (:,1) double
        low_hz  (1,1) double
        high_hz (1,1) double
    end

    mask = freq_hz >= low_hz & freq_hz <= high_hz;

    if sum(mask) < 2
        power_db = NaN;
        std_db   = NaN;
        return;
    end

    band_power = trapz(freq_hz(mask), psd(mask));
    power_db = 10 * log10(band_power + eps);
    std_db = std(psd_db(mask));
end
