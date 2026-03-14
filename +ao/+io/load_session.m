function session = load_session(case_dir)
%LOAD_SESSION  Load all MPX files in a case directory into a session struct.
%
%   session = ao.io.load_session('/path/to/case')
%
%   Discovers all .mpx files, parses filenames for depth/trajectory info,
%   reads each file, and groups them by trajectory.  Multiple files at the
%   same depth are concatenated.
%
%   Returns:
%       session.case_dir      - path to the case directory
%       session.segments()    - struct array, one per unique depth, with:
%           .filepath         - source .mpx file(s)
%           .hemisphere       - 'L' or 'R'
%           .trajectory       - trajectory number
%           .depth_mm         - depth in mm
%           .streams          - struct of continuous streams (from read_mpx)
%           .spikes           - struct of spike data (from read_mpx)

    arguments
        case_dir (1,:) char
    end

    if ~isfolder(case_dir)
        error('ao:io:load_session', 'Directory not found: %s', case_dir);
    end

    % --- Discover MPX files ---
    mpx_files = dir(fullfile(case_dir, '*.mpx'));
    if isempty(mpx_files)
        error('ao:io:load_session', 'No .mpx files found in: %s', case_dir);
    end

    % --- Parse filenames and sort ---
    n = numel(mpx_files);
    file_info = struct('filepath', cell(1,n), 'info', cell(1,n));
    valid_mask = false(1, n);

    for k = 1:n
        fp = fullfile(case_dir, mpx_files(k).name);
        info = ao.io.parse_filename(mpx_files(k).name);
        file_info(k).filepath = fp;
        file_info(k).info = info;
        valid_mask(k) = info.valid;
    end

    file_info = file_info(valid_mask);
    if isempty(file_info)
        error('ao:io:load_session', 'No valid MPX filenames found in: %s', case_dir);
    end

    % Sort by depth (ascending)
    depths = arrayfun(@(f) f.info.depth_mm, file_info);
    fidxs  = arrayfun(@(f) f.info.file_index, file_info);
    [~, order] = sortrows([depths(:), fidxs(:)]);
    file_info = file_info(order);

    % --- Group by (hemisphere, trajectory, depth) ---
    % Build a group key for each file
    group_keys = cell(numel(file_info), 1);
    for k = 1:numel(file_info)
        fi = file_info(k).info;
        group_keys{k} = sprintf('%s_T%d_D%.4f', fi.hemisphere, fi.trajectory, fi.depth_mm);
    end
    [unique_keys, ~, group_idx] = unique(group_keys, 'stable');

    % --- Read and assemble segments ---
    segments = struct([]);
    fprintf('Loading %d MPX files across %d depths...\n', numel(file_info), numel(unique_keys));

    for g = 1:numel(unique_keys)
        members = find(group_idx == g);
        first_info = file_info(members(1)).info;

        % Read all files in this group
        all_mpx = cell(numel(members), 1);
        all_paths = cell(numel(members), 1);
        for m = 1:numel(members)
            fp = file_info(members(m)).filepath;
            fprintf('  [%d/%d] %s\n', members(m), numel(file_info), fp);
            all_mpx{m} = ao.io.read_mpx(fp);
            all_paths{m} = fp;
        end

        % If multiple files at same depth, merge streams
        if numel(all_mpx) == 1
            merged_streams = all_mpx{1}.streams;
            merged_spikes  = all_mpx{1}.spikes;
        else
            merged_streams = merge_streams(all_mpx);
            merged_spikes  = all_mpx{1}.spikes;  % take first file's spikes
        end

        seg.filepath    = all_paths;
        seg.hemisphere  = first_info.hemisphere;
        seg.trajectory  = first_info.trajectory;
        seg.depth_mm    = first_info.depth_mm;
        seg.streams     = merged_streams;
        seg.spikes      = merged_spikes;

        if isempty(segments)
            segments = seg;
        else
            segments(end+1) = seg; %#ok<AGROW>
        end
    end

    session.case_dir = case_dir;
    session.segments = segments;
    fprintf('Loaded %d segments.\n', numel(segments));
end

%% ========================================================================
function merged = merge_streams(mpx_list)
%MERGE_STREAMS  Concatenate continuous streams from multiple MPX files.
    merged = mpx_list{1}.streams;
    stream_names = fieldnames(merged);

    for s = 1:numel(stream_names)
        sname = stream_names{s};
        for k = 2:numel(mpx_list)
            if isfield(mpx_list{k}.streams, sname)
                other = mpx_list{k}.streams.(sname);
                merged.(sname).data = [merged.(sname).data; other.data];
            end
        end
    end
end
