function plot_lfp_heatmap(ax, depths_mm, freq_hz, psd_db_matrix, band_limits, cfg)
%PLOT_LFP_HEATMAP  LFP spectral power heatmap (frequency x depth).
%
%   ao.viz.plot_lfp_heatmap(ax, depths, freq, psd_db_matrix, band_limits, cfg)
%
%   Uses pcolor (not imagesc) so non-uniform depth spacing is rendered
%   correctly — matching the Python pcolormesh behavior.

    arguments
        ax
        depths_mm     (:,1) double
        freq_hz       (:,1) double
        psd_db_matrix (:,:) double
        band_limits   (1,2) double
        cfg           struct
    end

    N = numel(depths_mm);

    % Build depth bin edges for pcolor (N+1 edges for N rows)
    if N > 1
        % Use midpoints between consecutive depths as edges
        mids = (depths_mm(1:end-1) + depths_mm(2:end)) / 2;
        depth_edges = [depths_mm(1) - (mids(1) - depths_mm(1)); ...
                       mids; ...
                       depths_mm(end) + (depths_mm(end) - mids(end))];
    else
        depth_edges = [depths_mm(1) - 0.5; depths_mm(1) + 0.5];
    end

    % Build frequency bin edges (N_freq+1 edges)
    df = median(diff(freq_hz));
    freq_edges = [freq_hz - df/2; freq_hz(end) + df/2];

    % Expand PSD matrix to (N+1) x (N_freq+1) for pcolor
    % pcolor drops the last row and column, so we pad
    psd_padded = [psd_db_matrix, psd_db_matrix(:, end)];
    psd_padded = [psd_padded; psd_padded(end, :)];

    [F, D] = meshgrid(freq_edges, depth_edges);
    pcolor(ax, F, D, psd_padded);
    shading(ax, 'flat');
    colormap(ax, parula);

    % Band highlight overlay
    hold(ax, 'on');
    yl = [min(depth_edges), max(depth_edges)];
    patch(ax, [band_limits(1) band_limits(2) band_limits(2) band_limits(1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0], ...
        'FaceAlpha', 0.15, 'EdgeColor', [0.9 0.7 0], 'LineWidth', 1);
    hold(ax, 'off');

    % Colorbar as compact legend in NE corner
    cb = colorbar(ax, 'Location', 'eastoutside');
    cb.Label.String = 'PSD (dB)';
    cb.FontSize = 6;

    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'Spectral Power vs Depth', 'FontSize', 10);
    set(ax, 'YDir', 'normal', 'FontSize', 8);
    xlim(ax, [cfg.lfp.fmin_hz, cfg.lfp.fmax_hz]);
end
