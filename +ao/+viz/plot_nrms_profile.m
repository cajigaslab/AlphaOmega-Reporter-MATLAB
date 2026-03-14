function plot_nrms_profile(ax, depths_mm, nrms_values, impedance_flags)
%PLOT_NRMS_PROFILE  Normalized RMS vs depth with impedance markers.
%
%   ao.viz.plot_nrms_profile(ax, depths, nrms, impedance_flags)

    arguments
        ax
        depths_mm       (:,1) double
        nrms_values     (:,1) double
        impedance_flags (:,1) logical = false(size(depths_mm))
    end

    hold(ax, 'on');

    % Plot valid depths
    valid = ~impedance_flags & isfinite(nrms_values);
    barh(ax, depths_mm(valid), nrms_values(valid), 0.6, ...
        'FaceColor', [0.2 0.6 0.9], 'EdgeColor', 'none');

    % Mark impedance depths with red X
    imp_idx = find(impedance_flags);
    for k = 1:numel(imp_idx)
        plot(ax, 0.5, depths_mm(imp_idx(k)), 'rx', ...
            'MarkerSize', 8, 'LineWidth', 2);
    end

    % Reference line at NRMS = 1
    xline(ax, 1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

    xlabel(ax, 'NRMS');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'NRMS vs Depth', 'FontSize', 10);
    set(ax, 'YDir', 'reverse', 'FontSize', 8);
    xlim(ax, [0, max(nrms_values(valid)) * 1.2 + 0.1]);
    hold(ax, 'off');
end
