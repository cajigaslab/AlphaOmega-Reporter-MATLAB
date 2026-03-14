function plot_nrms_profile(ax, depths_mm, nrms_values, impedance_flags)
%PLOT_NRMS_PROFILE  Normalized RMS vs depth — matches Python viz exactly.
%
%   ao.viz.plot_nrms_profile(ax, depths, nrms, impedance_flags)

    arguments
        ax
        depths_mm       (:,1) double
        nrms_values     (:,1) double
        impedance_flags (:,1) logical = false(size(depths_mm))
    end

    hold(ax, 'on');
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.2);
    set(ax, 'YDir', 'normal', 'FontSize', 9, 'Box', 'off');
    ax.XAxis.FontSize = 9;
    ax.YAxis.FontSize = 9;

    % Python: #355070 = [0.208 0.314 0.439]
    nrms_color = [0.208 0.314 0.439];

    % Plot valid (non-impedance) depths as connected line with markers
    valid = ~impedance_flags & isfinite(nrms_values);
    if any(valid)
        plot(ax, nrms_values(valid), depths_mm(valid), '-o', ...
            'Color', nrms_color, ...
            'MarkerFaceColor', nrms_color, ...
            'MarkerSize', 4, 'LineWidth', 1.4, ...
            'DisplayName', 'NRMS');
    end

    % Baseline line at NRMS = 1.0 — Python: #888888, --, lw=0.8
    xline(ax, 1, '--', 'Color', [0.533 0.533 0.533], 'LineWidth', 0.8, ...
        'DisplayName', 'Baseline (1.0)');

    % Mark impedance depths with red X at x=0 — Python: color='red', s=60
    if any(impedance_flags)
        imp_depths = depths_mm(impedance_flags);
        scatter(ax, zeros(size(imp_depths)), imp_depths, 60, 'rx', ...
            'LineWidth', 2, 'DisplayName', 'Impedance');
    end

    xlabel(ax, 'NRMS');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'NRMS vs Depth', 'FontSize', 10);

    legend(ax, 'Location', 'best', 'FontSize', 7);

    hold(ax, 'off');
end
