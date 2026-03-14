function cfg = defaults()
%DEFAULTS  Return the default report configuration struct.
%
%   cfg = ao.config.defaults()
%
%   All fields can be overridden by the user before passing to ao_report().

    % --- Frequency bands ---
    cfg.bands.delta.low_hz     = 1;    cfg.bands.delta.high_hz    = 4;
    cfg.bands.theta.low_hz     = 4;    cfg.bands.theta.high_hz    = 8;
    cfg.bands.alpha.low_hz     = 8;    cfg.bands.alpha.high_hz    = 13;
    cfg.bands.beta.low_hz      = 13;   cfg.bands.beta.high_hz     = 30;
    cfg.bands.highbeta.low_hz  = 20;   cfg.bands.highbeta.high_hz = 30;
    cfg.bands.gamma.low_hz     = 30;   cfg.bands.gamma.high_hz    = 100;

    cfg.primary_band = 'beta';
    cfg.profile_bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

    % --- LFP / spectral ---
    cfg.lfp.spectrum_method       = 'multitaper';  % 'multitaper' or 'welch'
    cfg.lfp.fmin_hz               = 0;
    cfg.lfp.fmax_hz               = 100;
    cfg.lfp.welch_nperseg         = 512;
    cfg.lfp.welch_noverlap        = 256;
    cfg.lfp.multitaper_nw         = 4.0;
    cfg.lfp.noise_notch_hz        = 60;    % set to [] to disable
    cfg.lfp.noise_notch_width_hz  = 2.0;

    % --- Spike detection ---
    cfg.detect.threshold_mad   = 4.5;
    cfg.detect.peak_sign       = 'neg';   % 'neg', 'pos', 'both'
    cfg.detect.dead_time_ms    = 0.8;
    cfg.detect.filter_hp_hz    = 300;
    cfg.detect.filter_lp_hz    = 6000;
    cfg.detect.filter_order    = 3;

    % --- Impedance check detection ---
    cfg.impedance.rms_ratio_threshold       = 10.0;
    cfg.impedance.kurtosis_threshold        = 20.0;
    cfg.impedance.spectral_entropy_threshold = 0.4;
    cfg.impedance.min_criteria              = 2;

    % --- RMS outlier detection ---
    cfg.outlier.modified_z_threshold = 3.5;

    % --- Rendering ---
    cfg.render.dpi                = 150;
    cfg.render.max_depth_mm       = 12.0;   % set to Inf to show all
    cfg.render.mer_trace_window_s = 0.250;
    cfg.render.raster_window_s    = 1.0;
    cfg.render.page_width_in      = 11;
    cfg.render.page_height_in     = 8.5;

    % --- Stream selection hints ---
    cfg.stream.mer_hints = {'SPK', 'RAW', 'MER', 'MICRO'};
    cfg.stream.lfp_hints = {'LFP', 'MACRO'};
end
