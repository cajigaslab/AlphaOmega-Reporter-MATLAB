# AlphaOmega-Reporter-MATLAB

A pure MATLAB implementation of the [AlphaOmega Reporter](https://github.com/cajigaslab/AlphaOmega-Reporter) for intraoperative neurophysiology recordings. Reads native `.mpx` files directly — no Python, no external dependencies beyond the Signal Processing Toolbox.

## What It Does

`ao_report` takes a folder of `.mpx` files recorded by an AlphaOmega Neuro Omega system during deep brain stimulation (DBS) surgery and produces a trajectory-centric PDF report showing the electrophysiological profile at every recorded depth.

### Report Panels

| Panel | Description |
|-------|-------------|
| **NRMS vs Depth** | Normalized RMS profile. Impedance-check depths marked with red X. |
| **Filtered MER vs Depth** | Spike-band filtered traces with detection thresholds and spike markers. |
| **Spectral Power vs Depth** | LFP power spectral density heatmap (frequency x depth) via multitaper estimation. |
| **Band Power vs Depth** | Mean power in delta, theta, alpha, beta, and gamma bands with +/-1 SD error bars. |

### Automatic Quality Control

- **Impedance check detection** — Composite 2-of-3 scoring (RMS ratio > 10x, kurtosis > 20, spectral entropy < 0.4). Flagged depths excluded from all panels and normalization.
- **RMS outlier flagging** — Modified Z-score (threshold 3.5) on non-impedance depths.

## Requirements

- MATLAB R2020b or later
- Signal Processing Toolbox (for `dpss`, `pwelch`, `butter`, `filtfilt`, `iirnotch`)

## Quickstart

```matlab
% Generate a report from the included sample case
ao_report('examples/sample_case')

% Specify output path
ao_report('/path/to/case', 'out', 'my_report.pdf')

% Custom configuration
cfg = ao.config.defaults();
cfg.lfp.spectrum_method = 'welch';
cfg.detect.threshold_mad = 5.0;
cfg.render.max_depth_mm = 15;
ao_report('/path/to/case', 'config', cfg)
```

## Project Structure

```
ao_report.m                 % Main entry point
+ao/
  +io/
    read_mpx.m              % Native MPX binary parser
    parse_filename.m         % Depth/trajectory extraction from filenames
    load_session.m           % Directory scanner and session assembler
    resolve_streams.m        % MER/LFP channel identification
  +features/
    compute_psd.m            % Multitaper and Welch PSD estimation
    compute_features.m       % RMS, kurtosis, entropy, Hjorth parameters
    compute_bandpower.m      % Integrated band power with std
    apply_notch_filter.m     % 60 Hz notch filter
    bandpass_filter.m        % Butterworth bandpass
  +detection/
    detect_spikes.m          % Threshold-based spike detection
  +qc/
    detect_impedance.m       % Impedance check detection (2-of-3 scoring)
    flag_rms_outliers.m      % Modified Z-score outlier detection
  +viz/
    plot_nrms_profile.m      % NRMS bar chart
    plot_mer_traces.m        % Stacked MER trace display
    plot_lfp_heatmap.m       % Spectral power heatmap
    plot_bandpower_profile.m % Multi-band power with error bars
  +config/
    defaults.m               % Default configuration struct
examples/
  sample_case/               % Anonymized 16-depth GPi DBS case
```

## Native MPX Parser

The MPX binary parser (`ao.io.read_mpx`) reads AlphaOmega format version 4 files directly. It handles:

- File header with recording metadata (date, time, version)
- Continuous analog channel definitions and data blocks
- Segmented (spike) channel definitions and waveform data
- Multi-file recordings at the same depth (automatic concatenation)
- Raw int16 to microvolt scaling using stored gain and amplitude parameters

No Python, no `neo` library, no external converters needed.

## Configuration

All parameters are controlled through the config struct returned by `ao.config.defaults()`:

| Category | Key Parameters |
|----------|---------------|
| **Spectral** | `lfp.spectrum_method` ('multitaper'/'welch'), `lfp.multitaper_nw` (4.0), `lfp.noise_notch_hz` (60) |
| **Bands** | `bands.beta.low_hz` (13), `bands.beta.high_hz` (30), `profile_bands` ({'delta','theta','alpha','beta','gamma'}) |
| **Detection** | `detect.threshold_mad` (4.5), `detect.peak_sign` ('neg'), `detect.dead_time_ms` (0.8) |
| **Impedance** | `impedance.rms_ratio_threshold` (10), `impedance.kurtosis_threshold` (20), `impedance.spectral_entropy_threshold` (0.4) |
| **Rendering** | `render.max_depth_mm` (12), `render.mer_trace_window_s` (0.250) |

## Relationship to the Python Version

This is a complete, independent rewrite. The Python version ([AlphaOmega-Reporter](https://github.com/cajigaslab/AlphaOmega-Reporter)) offers additional features including spike sorting, boundary decoding (HMM), multi-page reports, HDF5 export, and a macOS drag-and-drop app. This MATLAB version focuses on the core single-page report with native MPX reading for environments where MATLAB is the primary tool.

## License

See [LICENSE](LICENSE).
