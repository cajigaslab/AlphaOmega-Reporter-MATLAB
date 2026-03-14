function [mer_name, lfp_name] = resolve_streams(streams, cfg)
%RESOLVE_STREAMS  Find MER and LFP stream names from hints.
%
%   [mer_name, lfp_name] = ao.io.resolve_streams(streams, cfg)
%
%   Searches stream field names against cfg.stream.mer_hints and
%   cfg.stream.lfp_hints.  Returns '' if not found.

    arguments
        streams struct
        cfg     struct
    end

    names = fieldnames(streams);
    mer_name = '';
    lfp_name = '';

    % Find MER stream: match hints, prefer highest fs
    best_mer_fs = 0;
    for k = 1:numel(names)
        nm_upper = upper(names{k});
        s = streams.(names{k});
        for h = 1:numel(cfg.stream.mer_hints)
            if contains(nm_upper, upper(cfg.stream.mer_hints{h}))
                if s.fs_hz > best_mer_fs
                    mer_name = names{k};
                    best_mer_fs = s.fs_hz;
                end
                break;
            end
        end
        % Fallback: high sample rate likely MER
        if isempty(mer_name) && s.fs_hz >= 20000
            mer_name = names{k};
            best_mer_fs = s.fs_hz;
        end
    end

    % Find LFP stream: match hints, prefer lowest fs
    best_lfp_fs = Inf;
    for k = 1:numel(names)
        nm_upper = upper(names{k});
        s = streams.(names{k});
        for h = 1:numel(cfg.stream.lfp_hints)
            if contains(nm_upper, upper(cfg.stream.lfp_hints{h}))
                if s.fs_hz < best_lfp_fs
                    lfp_name = names{k};
                    best_lfp_fs = s.fs_hz;
                end
                break;
            end
        end
    end

    % Fallback: if only one stream, use it for both
    if numel(names) == 1
        if isempty(mer_name)
            mer_name = names{1};
        end
        if isempty(lfp_name)
            lfp_name = names{1};
        end
    end
end
