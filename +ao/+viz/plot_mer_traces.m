function plot_mer_traces(ax, depths_mm, time_s_cell, filtered_cell, ...
    threshold_uv_cell, spike_times_cell, impedance_flags, cfg)
%PLOT_MER_TRACES  Stacked filtered MER traces — matches Python viz exactly.
%
%   ao.viz.plot_mer_traces(ax, depths, times, filtered, thresholds, ...
%       spike_times, impedance_flags, cfg)
%
%   X-axis in seconds (matching Python). Traces auto-scaled per depth.

    arguments
        ax
        depths_mm        (:,1) double
        time_s_cell      (:,1) cell
        filtered_cell    (:,1) cell
        threshold_uv_cell (:,1) cell
        spike_times_cell  (:,1) cell
        impedance_flags  (:,1) logical
        cfg              struct
    end

    hold(ax, 'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.2);
    set(ax, 'YDir', 'normal', 'FontSize', 9, 'Box', 'off');

    N = numel(depths_mm);

    % --- Compute display gain (matches Python _trace_display_span_mm) ---
    visible = sort(unique(depths_mm));
    if numel(visible) >= 2
        spacings = diff(visible);
        spacings = spacings(spacings > 0 & isfinite(spacings));
        if ~isempty(spacings)
            depth_range = visible(end) - visible(1);
            median_based = median(spacings) * 0.5;
            range_based  = depth_range / max(numel(visible), 1) * 0.6;
            display_gain = max(median_based, range_based);
            display_gain = max(0.08, min(0.55, display_gain));
        else
            display_gain = 0.35;
        end
    else
        display_gain = 0.35;
    end

    % Python colors
    trace_color     = [0.122 0.122 0.122];  % #1f1f1f
    threshold_color = [0.757 0.071 0.122];  % #c1121f
    spike_color     = [0.043 0.443 0.537];  % #0b7189

    max_time = 0;

    for k = 1:N
        if impedance_flags(k)
            continue;
        end

        t = time_s_cell{k};
        filt = filtered_cell{k};
        if isempty(t) || isempty(filt)
            continue;
        end

        % Clip to display window
        win = cfg.render.mer_trace_window_s;
        mask = t <= win;
        t = t(mask);
        filt = filt(mask);

        if isempty(t)
            continue;
        end

        % Center the trace (subtract median, matching Python)
        filt = filt - median(filt, 'omitnan');

        % Scale: 97.5th percentile (matching Python _trace_display_scale)
        scale_uv = prctile(abs(filt), 97.5);
        if ~isfinite(scale_uv) || scale_uv < 1
            scale_uv = 1;
        end

        % Also consider threshold as scale candidate (Python behavior)
        if ~isempty(threshold_uv_cell{k}) && isfinite(threshold_uv_cell{k})
            scale_uv = max(scale_uv, abs(threshold_uv_cell{k}));
        end

        % Plot trace offset by depth — time in seconds (Python uses seconds)
        y_plot = depths_mm(k) + (filt / max(scale_uv, eps)) * display_gain;
        plot(ax, t, y_plot, 'Color', trace_color, 'LineWidth', 0.8);
        max_time = max(max_time, t(end));

        % Threshold lines — Python: #c1121f, --, lw=0.8, alpha=0.75
        if ~isempty(threshold_uv_cell{k}) && isfinite(threshold_uv_cell{k})
            thr = threshold_uv_cell{k};
            offset = (thr / max(scale_uv, eps)) * display_gain;
            % Negative threshold line (default peak_sign = 'neg')
            y_thr = depths_mm(k) - offset;
            line_h = plot(ax, [t(1) t(end)], [y_thr y_thr], '--', ...
                'Color', threshold_color, 'LineWidth', 0.8);
            line_h.Color(4) = 0.75;  % alpha
        end

        % Spike markers — Python: #0b7189, s=9
        sp_times = spike_times_cell{k};
        if ~isempty(sp_times)
            sp_in_win = sp_times(sp_times >= t(1) & sp_times <= t(end));
            if ~isempty(sp_in_win)
                sp_idx = interp1(t, 1:numel(t), sp_in_win, 'nearest', 'extrap');
                sp_idx = max(1, min(round(sp_idx), numel(filt)));
                sp_y = depths_mm(k) + (filt(sp_idx) / max(scale_uv, eps)) * display_gain;
                scatter(ax, sp_in_win, sp_y, 9, spike_color, 'filled', ...
                    'MarkerEdgeColor', 'none');
            end
        end
    end

    % Python: x-axis in seconds, xlabel "Time (s)"
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'Filtered MER vs Depth', 'FontSize', 10);
    xlim(ax, [0, max(max_time, 0.05)]);

    hold(ax, 'off');
end
