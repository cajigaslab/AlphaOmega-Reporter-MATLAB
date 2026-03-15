function [out_pdf, results] = ao_report(case_dir, varargin)
%AO_REPORT  Generate an AlphaOmega MER trajectory report.
%
%   ao_report('/path/to/case')
%   ao_report('/path/to/case', 'out', 'report.pdf')
%   ao_report('/path/to/case', 'config', cfg)
%   [pdf, results] = ao_report(...)
%
%   Reads all .mpx files in case_dir, computes per-depth features
%   (RMS, PSD, band power, spike detection, impedance checks), and
%   produces a landscape PDF report with four synchronized panels.
%   Each trajectory (hemisphere + trajectory number) gets its own page.
%
%   The second output RESULTS is a struct array (one element per
%   trajectory) containing all computed data.  It is also saved as
%   a .mat file alongside the PDF.

    p = inputParser;
    addRequired(p, 'case_dir', @ischar);
    addParameter(p, 'out', '', @ischar);
    addParameter(p, 'config', ao.config.defaults(), @isstruct);
    parse(p, case_dir, varargin{:});

    cfg = p.Results.config;
    out_pdf = p.Results.out;
    if isempty(out_pdf)
        [~, folder_name] = fileparts(case_dir);
        out_pdf = fullfile(case_dir, [folder_name '_report.pdf']);
    end

    % =====================================================================
    % 1. LOAD SESSION
    % =====================================================================
    session = ao.io.load_session(case_dir);
    segments = session.segments;

    if numel(segments) == 0
        error('ao:report', 'No valid segments found in: %s', case_dir);
    end

    % =====================================================================
    % 2. GROUP SEGMENTS BY TRAJECTORY
    % =====================================================================
    traj_keys = cell(numel(segments), 1);
    for i = 1:numel(segments)
        traj_keys{i} = sprintf('%sT%d', segments(i).hemisphere, segments(i).trajectory);
    end
    [unique_trajs, ~, traj_idx] = unique(traj_keys, 'stable');
    n_trajs = numel(unique_trajs);

    fprintf('Found %d trajectory(s): %s\n', n_trajs, strjoin(unique_trajs, ', '));

    % =====================================================================
    % 3. PROCESS EACH TRAJECTORY AND GENERATE PAGES
    % =====================================================================
    [~, folder_name] = fileparts(case_dir);
    all_figs = gobjects(n_trajs, 1);

    % Pre-allocate results struct array
    results = struct([]);

    for ti = 1:n_trajs
        seg_idx = find(traj_idx == ti);
        traj_segs = segments(seg_idx);
        N = numel(traj_segs);

        fprintf('\n=== Trajectory %s (%d depths) ===\n', unique_trajs{ti}, N);

        % --- Per-depth analysis ---
        [depths_mm, mer_rms, kurtosis_vals, entropy_vals, ...
            psd_results, band_power, band_std, ...
            time_s_cell, filtered_cell, threshold_cell, spike_times_cell, ...
            traj_warnings] = analyze_depths(traj_segs, cfg);

        band_names = cfg.profile_bands;

        % --- Impedance detection ---
        [impedance_flags, traj_warnings] = detect_impedances(...
            depths_mm, mer_rms, kurtosis_vals, entropy_vals, cfg, traj_warnings);

        % --- Mask impedance depths ---
        for i = find(impedance_flags)'
            for b = 1:numel(band_names)
                band_power.(band_names{b})(i) = NaN;
                band_std.(band_names{b})(i)   = NaN;
            end
            psd_results{i} = [];
        end

        % --- RMS outlier detection ---
        [rms_outliers, traj_warnings] = detect_rms_outliers(...
            depths_mm, mer_rms, impedance_flags, cfg, traj_warnings);

        % --- Build LFP heatmap matrix ---
        [psd_db_matrix, ref_freq] = build_psd_matrix(psd_results, N);

        % --- Merge duplicate depths ---
        nrms = NaN(N, 1);
        [depths_mm, nrms, mer_rms, impedance_flags, ...
            band_power, band_std, psd_db_matrix, ...
            time_s_cell, filtered_cell, threshold_cell, spike_times_cell, ...
            rms_outliers, N] = merge_duplicate_depths(...
            depths_mm, nrms, mer_rms, impedance_flags, ...
            band_power, band_std, psd_db_matrix, ...
            time_s_cell, filtered_cell, threshold_cell, spike_times_cell, ...
            rms_outliers, band_names);

        % --- NRMS computation ---
        non_imp = find(~impedance_flags & isfinite(mer_rms));
        baseline_rms = median(mer_rms(non_imp), 'omitnan');
        if baseline_rms > 0
            nrms(non_imp) = mer_rms(non_imp) / baseline_rms;
        end

        % --- Apply depth cap ---
        visible = true(N, 1);
        if isfinite(cfg.render.max_depth_mm)
            visible = depths_mm <= cfg.render.max_depth_mm;
        end
        vis_idx = find(visible);

        % --- Store results for this trajectory ---
        r = struct();
        r.trajectory      = unique_trajs{ti};
        r.case_name       = folder_name;
        r.depths_mm       = depths_mm;
        r.n_depths        = N;

        % MER features
        r.mer_rms         = mer_rms;
        r.nrms            = nrms;
        r.kurtosis        = kurtosis_vals;
        r.spectral_entropy = entropy_vals;

        % Spike detection (per-depth cell arrays)
        r.time_s          = time_s_cell;
        r.filtered_mer    = filtered_cell;
        r.threshold_uv    = threshold_cell;
        r.spike_times_s   = spike_times_cell;

        % LFP spectral
        r.psd_freq_hz     = ref_freq;
        r.psd_db_matrix   = psd_db_matrix;
        r.band_power      = band_power;
        r.band_std        = band_std;
        r.band_names      = band_names;

        % QC flags
        r.impedance_flags = impedance_flags;
        r.rms_outliers    = rms_outliers;

        % Warnings
        r.warnings        = traj_warnings;

        % Config snapshot
        r.config          = cfg;

        if isempty(results)
            results = r;
        else
            results(end+1) = r; %#ok<AGROW>
        end

        % --- Generate figure ---
        fig = render_trajectory_page(depths_mm, nrms, impedance_flags, ...
            time_s_cell, filtered_cell, threshold_cell, spike_times_cell, ...
            psd_db_matrix, ref_freq, band_power, band_std, ...
            vis_idx, traj_warnings, folder_name, unique_trajs{ti}, N, cfg);

        all_figs(ti) = fig;
    end

    % =====================================================================
    % 4. EXPORT ALL PAGES TO SINGLE PDF
    % =====================================================================
    for ti = 1:n_trajs
        if ti == 1
            exportgraphics(all_figs(ti), out_pdf, 'ContentType', 'image', ...
                'Resolution', cfg.render.dpi);
        else
            exportgraphics(all_figs(ti), out_pdf, 'ContentType', 'image', ...
                'Resolution', cfg.render.dpi, 'Append', true);
        end
        close(all_figs(ti));
    end

    fprintf('\nReport saved: %s\n', out_pdf);

    % =====================================================================
    % 5. SAVE RESULTS TO .MAT
    % =====================================================================
    [pdf_dir, pdf_stem] = fileparts(out_pdf);
    out_mat = fullfile(pdf_dir, [pdf_stem '_results.mat']);
    save(out_mat, 'results', '-v7.3');
    fprintf('Results saved: %s\n', out_mat);
end

%% ========================================================================
%  HELPER: Per-depth analysis
%  ========================================================================
function [depths_mm, mer_rms, kurtosis_vals, entropy_vals, ...
    psd_results, band_power, band_std, ...
    time_s_cell, filtered_cell, threshold_cell, spike_times_cell, ...
    warnings] = analyze_depths(traj_segs, cfg)

    N = numel(traj_segs);
    band_names = cfg.profile_bands;

    depths_mm       = zeros(N, 1);
    mer_rms         = NaN(N, 1);
    kurtosis_vals   = NaN(N, 1);
    entropy_vals    = ones(N, 1);
    psd_results     = cell(N, 1);
    band_power      = struct();
    band_std        = struct();
    time_s_cell     = cell(N, 1);
    filtered_cell   = cell(N, 1);
    threshold_cell  = cell(N, 1);
    spike_times_cell = cell(N, 1);
    warnings = {};

    for b = 1:numel(band_names)
        band_power.(band_names{b}) = NaN(N, 1);
        band_std.(band_names{b})   = NaN(N, 1);
    end

    for i = 1:N
        seg = traj_segs(i);
        depths_mm(i) = seg.depth_mm;

        [mer_name, lfp_name] = ao.io.resolve_streams(seg.streams, cfg);

        % --- MER analysis ---
        if ~isempty(mer_name) && isfield(seg.streams, mer_name)
            mer_stream = seg.streams.(mer_name);
            mer_data = strip_nan_padding(mer_stream.data);
            mer_fs = mer_stream.fs_hz;

            if numel(mer_data) > 100
                feat = ao.features.compute_features(mer_data, mer_fs);
                mer_rms(i)       = feat.mer_rms;
                kurtosis_vals(i) = feat.kurtosis;
                entropy_vals(i)  = feat.spectral_entropy;

                det = ao.detection.detect_spikes(mer_data, mer_fs, cfg);
                time_vec = (0:numel(det.filtered)-1)' / mer_fs;
                time_s_cell{i}      = time_vec;
                filtered_cell{i}    = det.filtered;
                threshold_cell{i}   = det.threshold_uv;
                spike_times_cell{i} = det.spike_times_s;
            end
        end

        % --- LFP analysis ---
        if ~isempty(lfp_name) && isfield(seg.streams, lfp_name)
            lfp_stream = seg.streams.(lfp_name);
            lfp_data = strip_nan_padding(lfp_stream.data);
            lfp_fs = lfp_stream.fs_hz;

            if numel(lfp_data) > 128
                if ~isempty(cfg.lfp.noise_notch_hz)
                    lfp_data = ao.features.apply_notch_filter(lfp_data, lfp_fs, ...
                        cfg.lfp.noise_notch_hz, cfg.lfp.noise_notch_width_hz);
                end

                psd_res = ao.features.compute_psd(lfp_data, lfp_fs, cfg);
                psd_results{i} = psd_res;

                if ~isempty(psd_res.freq_hz)
                    for b = 1:numel(band_names)
                        bn = band_names{b};
                        lo = cfg.bands.(bn).low_hz;
                        hi = cfg.bands.(bn).high_hz;
                        [pw, sd] = ao.features.compute_bandpower(...
                            psd_res.freq_hz, psd_res.psd, psd_res.psd_db, lo, hi);
                        band_power.(bn)(i) = pw;
                        band_std.(bn)(i)   = sd;
                    end
                end
            else
                warnings{end+1} = sprintf(...
                    'LFP PSD failed at depth %.2f mm: signal too short', ...
                    seg.depth_mm); %#ok<AGROW>
            end
        end

        fprintf('  [%d/%d] depth=%.2f mm  rms=%.1f\n', i, N, seg.depth_mm, mer_rms(i));
    end
end

%% ========================================================================
%  HELPER: Impedance detection
%  ========================================================================
function [impedance_flags, warnings] = detect_impedances(...
    depths_mm, mer_rms, kurtosis_vals, entropy_vals, cfg, warnings)

    N = numel(depths_mm);
    impedance_flags = false(N, 1);

    finite_rms = mer_rms(isfinite(mer_rms));
    median_rms = 0;
    if ~isempty(finite_rms)
        median_rms = median(finite_rms);
    end

    imp_strs = {};
    for i = 1:N
        if isfinite(mer_rms(i))
            impedance_flags(i) = ao.qc.detect_impedance(...
                mer_rms(i), kurtosis_vals(i), entropy_vals(i), median_rms, cfg);
            if impedance_flags(i)
                imp_strs{end+1} = sprintf('%.2f', depths_mm(i)); %#ok<AGROW>
            end
        end
    end
    if ~isempty(imp_strs)
        warnings{end+1} = ['Impedance check detected at: ' strjoin(imp_strs, ', ')];
    end
end

%% ========================================================================
%  HELPER: RMS outlier detection
%  ========================================================================
function [rms_outliers, warnings] = detect_rms_outliers(...
    depths_mm, mer_rms, impedance_flags, cfg, warnings)

    N = numel(depths_mm);
    rms_outliers = false(N, 1);
    non_imp = find(~impedance_flags & isfinite(mer_rms));

    if numel(non_imp) > 3
        [~, ol] = ao.qc.flag_rms_outliers(mer_rms(non_imp), cfg.outlier.modified_z_threshold);
        rms_outliers(non_imp) = ol;
        outlier_depths = depths_mm(non_imp(ol));
        if ~isempty(outlier_depths)
            strs = arrayfun(@(d) sprintf('%.2f', d), outlier_depths, 'UniformOutput', false);
            warnings{end+1} = ['RMS outliers at: ' strjoin(strs, ', ')];
        end
    end
end

%% ========================================================================
%  HELPER: Build PSD heatmap matrix
%  ========================================================================
function [psd_db_matrix, ref_freq] = build_psd_matrix(psd_results, N)
    ref_freq = [];
    for i = 1:N
        if ~isempty(psd_results{i}) && ~isempty(psd_results{i}.freq_hz)
            ref_freq = psd_results{i}.freq_hz;
            break;
        end
    end

    psd_db_matrix = NaN(N, numel(ref_freq));
    if ~isempty(ref_freq)
        for i = 1:N
            if ~isempty(psd_results{i}) && ~isempty(psd_results{i}.freq_hz)
                psd_db_matrix(i, :) = interp1(psd_results{i}.freq_hz, ...
                    psd_results{i}.psd_db, ref_freq, 'linear', NaN);
            end
        end
    end
end

%% ========================================================================
%  HELPER: Render one trajectory page
%  ========================================================================
function fig = render_trajectory_page(depths_mm, nrms, impedance_flags, ...
    time_s_cell, filtered_cell, threshold_cell, spike_times_cell, ...
    psd_db_matrix, ref_freq, band_power, band_std, ...
    vis_idx, warnings, folder_name, traj_label, N, cfg)

    band_names = cfg.profile_bands;

    fprintf('Generating figure for %s...\n', traj_label);

    fig = figure('Visible', 'off', 'PaperPositionMode', 'auto', ...
        'Units', 'inches', 'Position', [0 0 cfg.render.page_width_in cfg.render.page_height_in]);
    set(fig, 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', 9);

    % --- Layout with Python column width ratios [0.35, 1.15, 1.00, 0.85] ---
    % Python: constrained_layout with gridspec width_ratios
    ratios = [0.35, 1.15, 1.00, 0.85];
    total  = sum(ratios);

    % Leave room for: left margin, inter-panel gaps, right margin, footer
    % Shrink panel 3 ratio to account for colorbar taking ~5% of its width
    margin_l = 0.06;  margin_r = 0.03;  gap = 0.05;
    footer_h = 0.08;  top_h = 0.08;
    body_w = 1 - margin_l - margin_r - gap * (numel(ratios) - 1);
    body_b = footer_h;
    body_t = 1 - top_h;
    body_h = body_t - body_b;

    col_w = ratios / total * body_w;
    ax_positions = zeros(4, 4);  % [left bottom width height]
    x = margin_l;
    for c = 1:4
        ax_positions(c, :) = [x, body_b, col_w(c), body_h];
        x = x + col_w(c) + gap;
    end

    % --- Panel 1: NRMS ---
    ax1 = axes(fig, 'Position', ax_positions(1, :));
    ao.viz.plot_nrms_profile(ax1, depths_mm(vis_idx), nrms(vis_idx), ...
        impedance_flags(vis_idx));

    % --- Panel 2: Filtered MER traces ---
    ax2 = axes(fig, 'Position', ax_positions(2, :));
    ao.viz.plot_mer_traces(ax2, depths_mm(vis_idx), ...
        time_s_cell(vis_idx), filtered_cell(vis_idx), ...
        threshold_cell(vis_idx), spike_times_cell(vis_idx), ...
        impedance_flags(vis_idx), cfg);

    % --- Panel 3: LFP spectral heatmap ---
    ax3 = axes(fig, 'Position', ax_positions(3, :));
    if ~isempty(ref_freq)
        primary_band = cfg.primary_band;
        bl = [cfg.bands.(primary_band).low_hz, cfg.bands.(primary_band).high_hz];
        ao.viz.plot_lfp_heatmap(ax3, depths_mm(vis_idx), ref_freq, ...
            psd_db_matrix(vis_idx, :), bl, cfg);
    else
        title(ax3, 'No LFP data', 'FontSize', 10);
    end

    % --- Panel 4: Band power profiles ---
    ax4 = axes(fig, 'Position', ax_positions(4, :));

    % Python _default_band_colors(): exact hex-to-RGB conversions
    %   delta=#3366cc  theta=#2a9d8f  alpha=#2b9348
    %   beta=#f77f00   gamma=#d62828  highbeta=#7b2cbf
    band_colors = struct( ...
        'delta',    [0.200 0.400 0.800], ...
        'theta',    [0.165 0.616 0.561], ...
        'alpha',    [0.169 0.576 0.282], ...
        'beta',     [0.969 0.498 0.000], ...
        'gamma',    [0.839 0.157 0.157], ...
        'highbeta', [0.482 0.173 0.749]);

    bd_arr = struct([]);
    for b = 1:numel(band_names)
        bn = band_names{b};
        bd_arr(b).name = bn;
        bd_arr(b).power_db = band_power.(bn)(vis_idx);
        bd_arr(b).std_db   = band_std.(bn)(vis_idx);
        if isfield(band_colors, bn)
            bd_arr(b).color = band_colors.(bn);
        else
            bd_arr(b).color = [0.5 0.5 0.5];
        end
    end
    ao.viz.plot_bandpower_profile(ax4, depths_mm(vis_idx), bd_arr, cfg);

    % --- Synchronize y-axes (matching Python _depth_limits, pad=0.35) ---
    all_visible_depths = depths_mm(vis_idx);
    pad_mm = 0.35;
    yl = [min(all_visible_depths) - pad_mm, max(all_visible_depths) + pad_mm];
    ylim(ax1, yl); ylim(ax2, yl); ylim(ax3, yl); ylim(ax4, yl);

    % Remove y-tick labels from panels 2-4 (Python: sharey with label only on first)
    set(ax2, 'YTickLabel', []);
    set(ax3, 'YTickLabel', []);
    set(ax4, 'YTickLabel', []);
    ylabel(ax2, ''); ylabel(ax3, ''); ylabel(ax4, '');

    % Add shared y-label on figure left edge (Python: fig.supylabel)
    annotation(fig, 'textbox', [0.0 0.15 0.03 0.7], ...
        'String', 'Depth (mm)', 'FontSize', 9, ...
        'EdgeColor', 'none', 'Rotation', 90, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off');

    % --- Footer: warnings (Python wraps long text across multiple lines) ---
    if ~isempty(warnings)
        % Join with newlines instead of pipe for better wrapping
        warning_text = strjoin(warnings, newline);
        annotation(fig, 'textbox', [0.04 0.0 0.94 footer_h - 0.01], ...
            'String', warning_text, 'FontSize', 6, ...
            'EdgeColor', 'none', 'Color', [0.6 0 0], ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none', 'FitBoxToText', 'off');
    end

    % --- Title (Python: fig.suptitle, fontsize=11) ---
    title_str = sprintf('%s  |  %s  (%d depths)', folder_name, traj_label, N);
    sgtitle(fig, title_str, 'FontSize', 11, 'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

%% ========================================================================
function values = strip_nan_padding(data)
%STRIP_NAN_PADDING  Remove leading/trailing NaN from a data vector.
    data = data(:);
    first_valid = find(isfinite(data), 1, 'first');
    last_valid  = find(isfinite(data), 1, 'last');
    if isempty(first_valid)
        values = data;
    else
        values = data(first_valid:last_valid);
    end
end

%% ========================================================================
function [depths_out, nrms_out, rms_out, imp_out, ...
    bp_out, bs_out, psd_out, ...
    t_out, f_out, th_out, sp_out, ...
    ol_out, N_out] = merge_duplicate_depths(...
    depths, nrms, rms, imp_flags, ...
    bp, bs, psd_db, ...
    t_cell, f_cell, th_cell, sp_cell, ...
    ol_flags, band_names)
%MERGE_DUPLICATE_DEPTHS  Average features across same-depth segments.

    [u_depths, ~, ic] = unique(depths, 'stable');
    M = numel(u_depths);

    if M == numel(depths)
        depths_out = depths; nrms_out = nrms; rms_out = rms;
        imp_out = imp_flags; bp_out = bp; bs_out = bs;
        psd_out = psd_db; t_out = t_cell; f_out = f_cell;
        th_out = th_cell; sp_out = sp_cell; ol_out = ol_flags;
        N_out = M;
        return;
    end

    fprintf('  Merging %d segments -> %d unique depths\n', numel(depths), M);

    depths_out = u_depths;
    nrms_out   = NaN(M, 1);
    rms_out    = NaN(M, 1);
    imp_out    = false(M, 1);
    ol_out     = false(M, 1);

    nf = size(psd_db, 2);
    psd_out = NaN(M, nf);
    t_out  = cell(M, 1);
    f_out  = cell(M, 1);
    th_out = cell(M, 1);
    sp_out = cell(M, 1);

    bp_out = struct();
    bs_out = struct();
    for b = 1:numel(band_names)
        bp_out.(band_names{b}) = NaN(M, 1);
        bs_out.(band_names{b}) = NaN(M, 1);
    end

    for j = 1:M
        members = find(ic == j);
        rms_out(j)  = mean(rms(members), 'omitnan');
        nrms_out(j) = mean(nrms(members), 'omitnan');
        imp_out(j) = any(imp_flags(members));
        ol_out(j)  = any(ol_flags(members));

        for b = 1:numel(band_names)
            bn = band_names{b};
            bp_out.(bn)(j) = mean(bp.(bn)(members), 'omitnan');
            bs_out.(bn)(j) = mean(bs.(bn)(members), 'omitnan');
        end

        if nf > 0
            psd_out(j, :) = mean(psd_db(members, :), 1, 'omitnan');
        end

        best = members(1);
        best_len = 0;
        for k = 1:numel(members)
            m = members(k);
            if ~isempty(f_cell{m}) && numel(f_cell{m}) > best_len
                best_len = numel(f_cell{m});
                best = m;
            end
        end
        t_out{j}  = t_cell{best};
        f_out{j}  = f_cell{best};
        th_out{j} = th_cell{best};
        sp_out{j} = sp_cell{best};
    end

    N_out = M;
end
