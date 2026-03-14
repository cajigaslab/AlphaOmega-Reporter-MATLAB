function info = parse_filename(filename)
%PARSE_FILENAME  Extract hemisphere, trajectory, depth from an MPX filename.
%
%   info = ao.io.parse_filename('LT1D-3.683F0001.mpx')
%
%   Returns a struct with fields:
%       .hemisphere   - 'L', 'R', or '' (unknown)
%       .trajectory   - trajectory number (double)
%       .depth_mm     - depth in mm (double, signed)
%       .file_index   - file index within that depth (double)
%       .valid        - true if parsing succeeded
%
%   Pattern: [i]?[L|R]T<traj>D<depth>F<index>.mpx

    arguments
        filename (1,:) char
    end

    [~, stem, ~] = fileparts(filename);

    % Regex: optional leading 'i', optional hemisphere, T<n>D<depth>F<index>
    pat = '(?i)(?:i)?(?:(?<hemi>[RL]))?T(?<traj>\d+)D(?<depth>[+-]?\d+(?:\.\d+)?)F(?<fidx>\d+)';
    tok = regexp(stem, pat, 'names');

    if isempty(tok)
        info.hemisphere = '';
        info.trajectory = NaN;
        info.depth_mm   = NaN;
        info.file_index = NaN;
        info.valid      = false;
        return;
    end

    info.hemisphere = upper(tok.hemi);
    info.trajectory = str2double(tok.traj);
    info.depth_mm   = str2double(tok.depth);
    info.file_index = str2double(tok.fidx);
    info.valid      = true;
end
