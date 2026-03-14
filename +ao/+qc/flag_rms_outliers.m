function [zscores, is_outlier] = flag_rms_outliers(rms_values, threshold)
%FLAG_RMS_OUTLIERS  Modified Z-score outlier detection on RMS values.
%
%   [zscores, is_outlier] = ao.qc.flag_rms_outliers(rms_values, 3.5)
%
%   Modified Z-score: z* = 0.6745 × (x - median) / MAD
%   Outlier if |z*| > threshold (default 3.5).

    arguments
        rms_values (:,1) double
        threshold  (1,1) double = 3.5
    end

    med = median(rms_values, 'omitnan');
    mad_val = median(abs(rms_values - med), 'omitnan');

    if mad_val < eps
        zscores = zeros(size(rms_values));
    else
        zscores = 0.6745 * (rms_values - med) / mad_val;
    end

    is_outlier = abs(zscores) > threshold;
end
