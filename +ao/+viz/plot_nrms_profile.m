function plot_nrms_profile(ax, depths_mm, nrms_values, impedance_flags)
%PLOT_NRMS_PROFILE  Normalized RMS vs depth as a connected line plot.
%
%   ao.viz.plot_nrms_profile(ax, depths, nrms, impedance_flags)
%
%   Matches the Python report style: connected line with dot markers,
%   dashed baseline at NRMS=1, red X for impedance depths.

    arguments
        ax
        depths_mm       (:,1) double
        nrms_values     (:,1) double
        impedance_flags (:,1) logical = false(size(depths_mm))
    end

    hold(ax, 'on');

    % Reference line at NRMS = 1
    xline(ax, 1, '--', 'Color', [0.7 0.5 0.5], 'LineWidth', 0.8);

    % Plot valid (non-impedance) depths as connected line
    valid = ~impedance_flags & isfinite(nrms_values);
    if any(valid)
        plot(ax, nrms_values(valid), depths_mm(valid), '-o', ...
            'Color', [0.17 0.40 0.58], ...
            'MarkerFaceColor', [0.17 0.40 0.58], ...
            'MarkerSize', 3, 'LineWidth', 1.2);
    end

    % Mark impedance depths with red X
    imp_idx = find(impedance_flags);
    for k = 1:numel(imp_idx)
        plot(ax, 0.2, depths_mm(imp_idx(k)), 'rx', ...
            'MarkerSize', 8, 'LineWidth', 2);
    end

    xlabel(ax, 'NRMS');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'NRMS vs Depth', 'FontSize', 10);
    set(ax, 'YDir', 'normal', 'FontSize', 8);

    if any(valid)
        max_val = max(nrms_values(valid));
        xlim(ax, [0, max(max_val * 1.2, 1.5)]);
    else
        xlim(ax, [0, 2]);
    end

    % Build legend dynamically based on what was plotted
    leg_labels = {'Baseline (1.0)'};
    if any(valid)
        leg_labels{end+1} = 'NRMS';
    end
    if any(impedance_flags)
        leg_labels{end+1} = 'Impedance';
    end
    legend(ax, leg_labels, 'FontSize', 6, 'Location', 'northeast');

    hold(ax, 'off');
end
