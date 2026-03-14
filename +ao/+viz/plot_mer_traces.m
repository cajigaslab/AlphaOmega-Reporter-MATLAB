function plot_mer_traces(ax, depths_mm, time_s_cell, filtered_cell, ...
    threshold_uv_cell, spike_times_cell, impedance_flags, cfg)
%PLOT_MER_TRACES  Stacked filtered MER traces with spike markers.
%
%   ao.viz.plot_mer_traces(ax, depths, times, filtered, thresholds, ...
%       spike_times, impedance_flags, cfg)
%
%   Each depth's trace is offset vertically by its depth value and
%   amplitude-scaled to fit within the inter-depth spacing.

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
    N = numel(depths_mm);

    % Compute display gain from median depth spacing
    if N > 1
        spacings = abs(diff(sort(depths_mm)));
        median_spacing = median(spacings(spacings > 0));
        display_gain = max(median_spacing * 0.5, ...
            (max(depths_mm) - min(depths_mm)) / N * 0.6);
        display_gain = max(0.08, min(0.55, display_gain));
    else
        display_gain = 0.4;
    end

    for k = 1:N
        if impedance_flags(k)
            continue;  % skip impedance depths
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

        % Scale: auto based on 97.5th percentile
        scale_uv = prctile(abs(filt), 97.5);
        if scale_uv < 1
            scale_uv = 1;
        end

        % Plot trace offset by depth
        y_plot = depths_mm(k) + (filt / scale_uv) * display_gain;
        plot(ax, t * 1000, y_plot, 'Color', [0.15 0.15 0.15], 'LineWidth', 0.4);

        % Threshold lines
        if ~isempty(threshold_uv_cell{k})
            thr = threshold_uv_cell{k};
            y_thr_neg = depths_mm(k) + (-thr / scale_uv) * display_gain;
            yline_at = y_thr_neg;
            plot(ax, [t(1) t(end)] * 1000, [yline_at yline_at], ...
                '--', 'Color', [0.9 0.2 0.2], 'LineWidth', 0.5);
        end

        % Spike markers
        sp_times = spike_times_cell{k};
        if ~isempty(sp_times)
            sp_in_win = sp_times(sp_times <= win);
            if ~isempty(sp_in_win)
                % Interpolate filtered signal at spike times for y position
                sp_vals = interp1(t, filt, sp_in_win, 'nearest', 'extrap');
                sp_y = depths_mm(k) + (sp_vals / scale_uv) * display_gain;
                plot(ax, sp_in_win * 1000, sp_y, '.', ...
                    'Color', [0 0.8 0.8], 'MarkerSize', 4);
            end
        end
    end

    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'Filtered MER vs Depth', 'FontSize', 10);
    set(ax, 'YDir', 'normal', 'FontSize', 8);
    hold(ax, 'off');
end
