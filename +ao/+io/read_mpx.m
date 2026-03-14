function session = read_mpx(filepath)
%READ_MPX  Read an AlphaOmega .mpx file into a MATLAB struct.
%
%   session = ao.io.read_mpx(filepath)
%
%   Parses the binary MPX format (version 4) natively — no Python or neo
%   dependency.  Returns a struct with fields:
%
%       .filepath       - source file path
%       .header         - file header metadata (date/time, version)
%       .channels       - struct array of channel definitions
%       .streams        - struct mapping channel names to continuous data
%       .spikes         - struct mapping channel names to spike/segment data
%       .digital        - struct mapping channel names to digital event data
%
%   Each stream entry contains:
%       .data       - double matrix (N_samples × N_channels), in µV
%       .fs_hz      - sampling rate in Hz
%       .t0_sample  - first sample number (hardware counter)
%       .name       - channel name string
%
%   Each spike entry contains:
%       .times_s    - spike times in seconds (N_spikes × 1)
%       .waveforms  - waveform matrix (N_spikes × N_samples), in µV
%       .unit_ids   - unit labels (N_spikes × 1)
%       .fs_hz      - waveform sampling rate in Hz
%       .pre_ms     - pre-trigger time in ms
%       .post_ms    - post-trigger time in ms
%       .name       - channel name string

    arguments
        filepath (1,:) char
    end

    fid = fopen(filepath, 'r', 'l');  % little-endian
    if fid == -1
        error('ao:io:read_mpx', 'Cannot open file: %s', filepath);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    % --- Read file header (first block, type 'h') ---
    hdr = read_header_block(fid);

    % --- Parse all blocks ---
    chan_defs = struct([]);   % channel definitions
    data_blocks = {};        % raw data block list
    digital_events = {};     % digital event list

    while true
        % Read block header: 2-byte length + 1-byte type
        raw = fread(fid, 3, '*uint8');
        if numel(raw) < 3
            break;  % EOF
        end
        block_len = typecast(raw(1:2), 'uint16');
        block_type = char(raw(3));

        if block_len == 65535  % 0xFFFF = EOF sentinel
            break;
        end

        % Read remaining bytes for this block (length includes the 3-byte header)
        remaining = int32(block_len) - 3;
        if remaining < 0
            break;
        end
        if remaining > 0
            block_data = fread(fid, remaining, '*uint8');
            if numel(block_data) < remaining
                break;  % truncated file
            end
        else
            block_data = uint8([]);
        end

        switch block_type
            case '2'
                def = parse_channel_def(block_data, block_len);
                if ~isempty(def)
                    if isempty(chan_defs)
                        chan_defs = def;
                    else
                        chan_defs(end+1) = def; %#ok<AGROW>
                    end
                end

            case '5'
                db = parse_data_block(block_data);
                if ~isempty(db)
                    data_blocks{end+1} = db; %#ok<AGROW>
                end

            case 'S'
                % Stream definition — skip for now
            case 'b'
                % Port definition — skip for now
            case 'E'
                % Event data — skip for now
            otherwise
                % Unknown block type — skip
        end
    end

    % --- Build channel lookup ---
    chan_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:numel(chan_defs)
        chan_map(int32(chan_defs(k).channel_number)) = chan_defs(k);
    end

    % --- Assemble continuous streams and spike data ---
    streams = struct();
    spikes  = struct();

    % Group data blocks by channel number
    block_groups = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:numel(data_blocks)
        db = data_blocks{k};
        ch = int32(db.channel_number);
        if block_groups.isKey(ch)
            block_groups(ch) = [block_groups(ch), {db}];
        else
            block_groups(ch) = {db};
        end
    end

    ch_nums = block_groups.keys();
    for k = 1:numel(ch_nums)
        ch_num = ch_nums{k};
        blocks = block_groups(ch_num);

        if ~chan_map.isKey(ch_num)
            continue;  % no definition for this channel
        end
        def = chan_map(ch_num);

        if def.mode == 0
            % --- Continuous channel ---
            s = assemble_continuous(blocks, def);
            fname = sanitize_fieldname(def.name);
            streams.(fname) = s;
        else
            % --- Segmented (spike) channel ---
            s = assemble_spikes(blocks, def);
            fname = sanitize_fieldname(def.name);
            spikes.(fname) = s;
        end
    end

    % --- Build output ---
    session.filepath = filepath;
    session.header   = hdr;
    session.channels = chan_defs;
    session.streams  = streams;
    session.spikes   = spikes;
end

%% ========================================================================
%  HEADER BLOCK
%  ========================================================================
function hdr = read_header_block(fid)
    raw = fread(fid, 3, '*uint8');
    block_len = typecast(raw(1:2), 'uint16');
    block_type = char(raw(3));
    assert(block_type == 'h', 'ao:io:read_mpx', ...
        'First block must be header (type h), got %s', block_type);

    remaining = int32(block_len) - 3;
    data = fread(fid, remaining, '*uint8');

    % Parse SDataHeader (56 bytes after alignment byte)
    % Byte 0: alignment (skip)
    idx = 2;  % 1-indexed, skip alignment byte at position 1
    hdr.next_datablock = typecast(data(idx:idx+3), 'uint32'); idx = idx + 4;
    hdr.version        = typecast(data(idx:idx+1), 'int16');  idx = idx + 2;
    hdr.hour           = data(idx); idx = idx + 1;
    hdr.minute         = data(idx); idx = idx + 1;
    hdr.second         = data(idx); idx = idx + 1;
    hdr.hsecond        = data(idx); idx = idx + 1;
    hdr.day            = data(idx); idx = idx + 1;
    hdr.month          = data(idx); idx = idx + 1;
    hdr.year           = typecast(data(idx:idx+1), 'uint16'); idx = idx + 2;
    hdr.dayofweek      = data(idx); idx = idx + 1;
    idx = idx + 1;  % padding
    hdr.minimum_time   = typecast(data(idx:idx+7), 'double'); idx = idx + 8;
    hdr.maximum_time   = typecast(data(idx:idx+7), 'double'); idx = idx + 8;
    hdr.erase_count    = typecast(data(idx:idx+3), 'uint32'); idx = idx + 4;
    hdr.map_version    = data(idx); idx = idx + 1;

    % Application name (10 bytes, null-terminated)
    app_raw = data(idx:idx+9); idx = idx + 10;
    null_pos = find(app_raw == 0, 1);
    if ~isempty(null_pos)
        hdr.application = char(app_raw(1:null_pos-1));
    else
        hdr.application = char(app_raw);
    end

    hdr.resource_version = data(idx:idx+3); idx = idx + 4; %#ok<NASGU>

    if hdr.map_version ~= 4
        warning('ao:io:read_mpx', ...
            'Expected MPX format version 4, got %d', hdr.map_version);
    end
end

%% ========================================================================
%  CHANNEL DEFINITION (Type 2)
%  ========================================================================
function def = parse_channel_def(data, block_len)
    def = [];
    if numel(data) < 12
        return;
    end

    idx = 1;
    % alignment byte
    idx = idx + 1;
    % next_datablock
    idx = idx + 4;  % skip (we don't use linked-list traversal)
    is_analog = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;
    is_input  = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;
    channel_number = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;

    % Skip alignment + RGB color (4 bytes)
    idx = idx + 4;

    if is_analog == 1 && is_input == 1
        % SDefAnalog (10 bytes)
        mode = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;
        amplitude = typecast(data(idx:idx+3), 'single'); idx = idx + 4;
        sample_rate_khz = typecast(data(idx:idx+3), 'single'); idx = idx + 4;

        % Fix amplitude if too small
        if amplitude <= 5
            amplitude = single(1250000 / 2^15);
        end

        spike_count = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;
        mode_spike  = typecast(data(idx:idx+1), 'uint16'); idx = idx + 2;

        if mode == 0
            % Continuous analog: SDefContinAnalog (6 bytes + name)
            if idx + 5 > numel(data)
                return;
            end
            duration = typecast(data(idx:idx+3), 'single'); idx = idx + 4;
            total_gain_100 = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;

            % Name: remaining bytes, null-terminated
            name = read_null_string(data, idx);

            def.channel_number = double(channel_number);
            def.is_analog = true;
            def.mode = 0;  % continuous
            def.amplitude = double(amplitude);
            def.sample_rate_hz = double(sample_rate_khz) * 1000;
            def.spike_count = double(spike_count);
            def.duration = double(duration);
            def.total_gain_100 = double(total_gain_100);
            def.name = name;
            def.pre_ms = 0;
            def.post_ms = 0;

        elseif mode == 1
            % Segmented analog: SDefLevelAnalog (14 bytes + name)
            if idx + 13 > numel(data)
                return;
            end
            pre_trig_ms  = typecast(data(idx:idx+3), 'single'); idx = idx + 4;
            post_trig_ms = typecast(data(idx:idx+3), 'single'); idx = idx + 4;
            idx = idx + 2;  % level_value
            idx = idx + 2;  % trg_mode
            idx = idx + 2;  % yes_rms
            total_gain_100 = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;

            name = read_null_string(data, idx);

            def.channel_number = double(channel_number);
            def.is_analog = true;
            def.mode = 1;  % segmented
            def.amplitude = double(amplitude);
            def.sample_rate_hz = double(sample_rate_khz) * 1000;
            def.spike_count = double(spike_count);
            def.duration = 0;
            def.total_gain_100 = double(total_gain_100);
            def.name = name;
            def.pre_ms = double(pre_trig_ms);
            def.post_ms = double(post_trig_ms);
        else
            return;
        end

    elseif is_analog == 0 && is_input == 1
        % Digital input channel
        if idx + 11 > numel(data)
            return;
        end
        sample_rate_khz = typecast(data(idx:idx+3), 'single'); idx = idx + 4;
        idx = idx + 2;  % save_trigger
        idx = idx + 4;  % duration
        idx = idx + 2;  % prev_status

        name = read_null_string(data, idx);

        def.channel_number = double(channel_number);
        def.is_analog = false;
        def.mode = -1;  % digital
        def.amplitude = 1;
        def.sample_rate_hz = double(sample_rate_khz) * 1000;
        def.spike_count = 0;
        def.duration = 0;
        def.total_gain_100 = 100;
        def.name = name;
        def.pre_ms = 0;
        def.post_ms = 0;
    else
        return;  % output channel or unknown
    end
end

%% ========================================================================
%  DATA BLOCK (Type 5)
%  ========================================================================
function db = parse_data_block(data)
    db = [];
    if numel(data) < 3
        return;
    end

    idx = 1;
    unit_number = data(idx); idx = idx + 1;
    channel_number = typecast(data(idx:idx+1), 'int16'); idx = idx + 2;

    % Remaining bytes: samples (int16) + trailing uint32 sample number
    sample_bytes = numel(data) - 3;  % after unit + channel_number

    if sample_bytes < 4
        return;  % need at least the trailing sample counter
    end

    % Last 4 bytes are the first_sample_number (uint32)
    first_sample = typecast(data(end-3:end), 'uint32');

    % Sample data: everything between channel_number and the trailing counter
    n_sample_bytes = sample_bytes - 4;
    if n_sample_bytes < 2
        db.channel_number = double(channel_number);
        db.unit_number = double(unit_number);
        db.samples = int16([]);
        db.first_sample = double(first_sample);
        return;
    end

    raw_samples = typecast(data(idx:idx+n_sample_bytes-1), 'int16');

    db.channel_number = double(channel_number);
    db.unit_number = double(unit_number);
    db.samples = raw_samples;
    db.first_sample = double(first_sample);
end

%% ========================================================================
%  ASSEMBLE CONTINUOUS STREAM
%  ========================================================================
function s = assemble_continuous(blocks, def)
    % Sort blocks by first_sample
    first_samples = cellfun(@(b) b.first_sample, blocks);
    [~, order] = sort(first_samples);
    blocks = blocks(order);

    % Compute scale factor: raw_int16 -> µV
    gain = def.total_gain_100 / 100;
    if gain == 0
        gain = 1;
    end
    scale = def.amplitude / gain;

    % Find total span
    all_samples_raw = cell(numel(blocks), 1);
    all_offsets = zeros(numel(blocks), 1);
    min_sample = blocks{1}.first_sample;

    for k = 1:numel(blocks)
        all_samples_raw{k} = blocks{k}.samples;
        all_offsets(k) = blocks{k}.first_sample - min_sample;
    end

    % Total length
    max_end = 0;
    for k = 1:numel(blocks)
        block_end = all_offsets(k) + numel(all_samples_raw{k});
        if block_end > max_end
            max_end = block_end;
        end
    end

    % Pre-allocate with NaN (gaps become NaN)
    assembled = nan(max_end, 1);
    for k = 1:numel(blocks)
        start_idx = all_offsets(k) + 1;
        n = numel(all_samples_raw{k});
        assembled(start_idx:start_idx+n-1) = double(all_samples_raw{k}) * scale;
    end

    s.data = assembled;
    s.fs_hz = def.sample_rate_hz;
    s.t0_sample = min_sample;
    s.name = def.name;
    s.channel_number = def.channel_number;
    s.amplitude = def.amplitude;
    s.total_gain_100 = def.total_gain_100;
end

%% ========================================================================
%  ASSEMBLE SPIKE / SEGMENTED DATA
%  ========================================================================
function s = assemble_spikes(blocks, def)
    fs_hz = def.sample_rate_hz;
    pre_samples  = round(def.pre_ms / 1000 * fs_hz);
    post_samples = round(def.post_ms / 1000 * fs_hz);
    wf_len = def.spike_count;
    if wf_len == 0
        wf_len = pre_samples + post_samples;
    end

    gain = def.total_gain_100 / 100;
    if gain == 0
        gain = 1;
    end
    scale = def.amplitude / gain;

    times   = zeros(numel(blocks), 1);
    units   = zeros(numel(blocks), 1);
    wfs     = nan(numel(blocks), max(wf_len, 1));

    valid = 0;
    for k = 1:numel(blocks)
        b = blocks{k};
        times_sample = b.first_sample;
        t_s = times_sample / fs_hz;

        n_samp = numel(b.samples);
        if n_samp == 0
            continue;
        end

        valid = valid + 1;
        times(valid) = t_s;
        units(valid) = b.unit_number;

        n_copy = min(n_samp, size(wfs, 2));
        wfs(valid, 1:n_copy) = double(b.samples(1:n_copy)) * scale;
    end

    s.times_s   = times(1:valid);
    s.waveforms = wfs(1:valid, :);
    s.unit_ids  = units(1:valid);
    s.fs_hz     = fs_hz;
    s.pre_ms    = def.pre_ms;
    s.post_ms   = def.post_ms;
    s.name      = def.name;
end

%% ========================================================================
%  HELPERS
%  ========================================================================
function str = read_null_string(data, start_idx)
    if start_idx > numel(data)
        str = '';
        return;
    end
    sub = data(start_idx:end);
    null_pos = find(sub == 0, 1);
    if isempty(null_pos)
        str = char(sub);
    else
        str = char(sub(1:null_pos-1));
    end
    str = strtrim(str);
end

function fname = sanitize_fieldname(name)
    % Convert channel name to a valid MATLAB struct field name
    fname = regexprep(name, '[^a-zA-Z0-9_]', '_');
    if isempty(fname) || ~isstrprop(fname(1), 'alpha')
        fname = ['ch_' fname];
    end
    fname = matlab.lang.makeValidName(fname);
end
