function plot_bandpower_profile(ax, depths_mm, band_data, cfg)
%PLOT_BANDPOWER_PROFILE  Multi-band power vs depth with ±1 SD error bars.
%
%   ao.viz.plot_bandpower_profile(ax, depths, band_data, cfg)
%
%   band_data: struct array with fields:
%       .name       - band name string (e.g. 'beta')
%       .power_db   - power values per depth (N×1, NaN for missing)
%       .std_db     - std values per depth (N×1, NaN for missing)
%       .color      - [R G B] color

    arguments
        ax
        depths_mm (:,1) double
        band_data (:,1) struct
        cfg       struct
    end

    hold(ax, 'on');
    legend_entries = {};
    legend_handles = [];

    for b = 1:numel(band_data)
        bd = band_data(b);
        valid = isfinite(bd.power_db);

        if ~any(valid)
            continue;
        end

        col = bd.color;
        light_col = col + (1 - col) * 0.4;  % lighten for error bars

        d = depths_mm(valid);
        p = bd.power_db(valid);
        s = bd.std_db(valid);
        s(~isfinite(s)) = 0;

        h = errorbar(ax, p, d, [], [], s, s, ...
            'Color', light_col, 'LineWidth', 0.8, 'CapSize', 2);
        h.LineStyle = 'none';

        % Plot mean line on top
        hline = plot(ax, p, d, '-o', 'Color', col, ...
            'LineWidth', 1.2, 'MarkerSize', 3, 'MarkerFaceColor', col);

        legend_handles(end+1) = hline; %#ok<AGROW>
        % Include frequency range in legend if available
        if isfield(cfg, 'bands') && isfield(cfg.bands, bd.name)
            lo = cfg.bands.(bd.name).low_hz;
            hi = cfg.bands.(bd.name).high_hz;
            legend_entries{end+1} = sprintf('%s (%g-%g Hz)', bd.name, lo, hi); %#ok<AGROW>
        else
            legend_entries{end+1} = bd.name; %#ok<AGROW>
        end
    end

    legend(ax, legend_handles, legend_entries, ...
        'FontSize', 7, 'Location', 'northeast');

    xlabel(ax, 'Power (dB)');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'Band Power vs Depth', 'FontSize', 10);
    set(ax, 'YDir', 'normal', 'FontSize', 8);

    % Add annotation
    text(ax, 0.02, 0.02, 'mean \pm 1 SD', ...
        'Units', 'normalized', 'FontSize', 6, ...
        'Color', [0.5 0.5 0.5]);

    hold(ax, 'off');
end
