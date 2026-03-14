function plot_lfp_heatmap(ax, depths_mm, freq_hz, psd_db_matrix, band_limits, cfg)
%PLOT_LFP_HEATMAP  LFP spectral power heatmap (frequency × depth).
%
%   ao.viz.plot_lfp_heatmap(ax, depths, freq, psd_db_matrix, band_limits, cfg)
%
%   psd_db_matrix: N_depths × N_freqs matrix (rows = depths, cols = freqs)
%   band_limits:   [low_hz, high_hz] for the primary band highlight

    arguments
        ax
        depths_mm     (:,1) double
        freq_hz       (:,1) double
        psd_db_matrix (:,:) double
        band_limits   (1,2) double
        cfg           struct
    end

    % Plot heatmap
    imagesc(ax, freq_hz, depths_mm, psd_db_matrix);
    colormap(ax, parula);

    % Add semi-transparent band highlight
    hold(ax, 'on');
    yl = [min(depths_mm) - 0.3, max(depths_mm) + 0.3];
    patch(ax, [band_limits(1) band_limits(2) band_limits(2) band_limits(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0], ...
        'FaceAlpha', 0.15, 'EdgeColor', [0.9 0.7 0], 'LineWidth', 1);
    hold(ax, 'off');

    % Colorbar as legend in NE corner
    cbar_pos = get(ax, 'Position');
    cb = colorbar(ax, 'Location', 'manual');
    cb_width  = cbar_pos(3) * 0.06;
    cb_height = cbar_pos(4) * 0.30;
    cb.Position = [cbar_pos(1) + cbar_pos(3) - cb_width - cbar_pos(3)*0.02, ...
                   cbar_pos(2) + cbar_pos(4) - cb_height - cbar_pos(4)*0.02, ...
                   cb_width, cb_height];
    cb.Label.String = 'dB';
    cb.FontSize = 6;
    cb.Color = [1 1 1];
    % Make semi-transparent background via box
    cb.Box = 'on';
    cb.TickDirection = 'out';

    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'Spectral Power vs Depth', 'FontSize', 10);
    set(ax, 'YDir', 'normal', 'FontSize', 8);
    xlim(ax, [cfg.lfp.fmin_hz, cfg.lfp.fmax_hz]);
end
