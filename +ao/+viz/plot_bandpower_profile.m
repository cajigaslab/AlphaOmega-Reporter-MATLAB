function plot_bandpower_profile(ax, depths_mm, band_data, cfg)
%PLOT_BANDPOWER_PROFILE  Multi-band power vs depth — matches Python viz exactly.
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
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.2);
    set(ax, 'YDir', 'normal', 'FontSize', 9, 'Box', 'off');

    legend_handles = [];
    legend_labels  = {};
    xmax = -Inf;
    xmin = Inf;

    for b = 1:numel(band_data)
        bd = band_data(b);
        valid = isfinite(bd.power_db);

        if ~any(valid)
            continue;
        end

        col = bd.color;
        % Lighten color for error bars (matching Python _lighten_color, factor=0.4)
        light_col = col + (1 - col) * 0.4;

        d = depths_mm(valid);
        p = bd.power_db(valid);
        s = bd.std_db(valid);
        s(~isfinite(s)) = 0;

        % Error bars — Python: elinewidth=0.8, capsize=0.0
        errorbar(ax, p, d, [], [], s, s, ...
            'Color', light_col, 'LineWidth', 0.8, 'CapSize', 0, ...
            'LineStyle', 'none');

        % Mean line on top — Python: fmt='o-', lw=1.2, markersize=3.5
        hline = plot(ax, p, d, '-o', 'Color', col, ...
            'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', col);

        legend_handles(end+1) = hline; %#ok<AGROW>

        % Include frequency range in legend
        if isfield(cfg, 'bands') && isfield(cfg.bands, bd.name)
            lo = cfg.bands.(bd.name).low_hz;
            hi = cfg.bands.(bd.name).high_hz;
            legend_labels{end+1} = sprintf('%s (%g-%g Hz)', bd.name, lo, hi); %#ok<AGROW>
        else
            legend_labels{end+1} = bd.name; %#ok<AGROW>
        end

        xmax = max(xmax, max(p));
        xmin = min(xmin, min(p));
    end

    if isempty(legend_handles)
        title(ax, 'Band Power vs Depth', 'FontSize', 10);
        hold(ax, 'off');
        return;
    end

    % X-axis padding — Python: span * 0.08, min 0.5
    span = xmax - xmin;
    if ~isfinite(span), span = 0; end
    pad = max(span * 0.08, 0.5);
    xlim(ax, [xmin - pad, xmax + pad]);

    % Legend — Python: loc='upper right', frameon=True, framealpha=0.9,
    %          handlelength=1.8, fontsize=7
    lg = legend(ax, legend_handles, legend_labels, ...
        'FontSize', 7, 'Location', 'northeast');
    lg.BoxFace.ColorType = 'truecoloralpha';
    lg.BoxFace.ColorData = uint8([255; 255; 255; 230]);  % framealpha ~0.9
    lg.EdgeColor = [0.816 0.816 0.816];  % #d0d0d0

    % Annotation — Python: "mean ± 1 SD", fontsize=7, lower-right
    text(ax, 0.98, 0.02, 'mean \pm 1 SD', ...
        'Units', 'normalized', 'FontSize', 7, ...
        'Color', [0.5 0.5 0.5], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    xlabel(ax, 'Bandpower (dB)');
    ylabel(ax, 'Depth (mm)');
    title(ax, 'Band Power vs Depth', 'FontSize', 10);
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 8;

    hold(ax, 'off');
end
