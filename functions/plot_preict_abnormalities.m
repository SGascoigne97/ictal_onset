function plot_preict_abnormalities(data_tbl, cell_imprint, opts)
% Plot imprints on EEG for determining performance of imprint algorithm
% input:
%   - data_tbl: full data table
%   - cell_imprint: table of imprints for individual subject
%   - optional inputs
%       - save_fig: logical determining if we want to store imprints
%       - save_fig_loc: minimum number of seizures to have been recorded per patient
%       - imprint_changes_label: string with brief description of how
%       imprint has been altered

    arguments
        data_tbl
        cell_imprint
        opts.save_fig = 1
        opts.save_fig_loc = "figures/"
        opts.imprint_changes_label = "original"
        opts.fig_visible (1,1) string {mustBeMember(opts.fig_visible, ["on", "off"])} = "on"
        opts.sz (1,1) double = NaN
        opts.channel_labels = NaN
    end

    % Set optional arguments
    save_fig = opts.save_fig;
    save_fig_loc = opts.save_fig_loc;
    imprint_changes_label = opts.imprint_changes_label; 
    fig_visible = opts.fig_visible;
    sz = opts.sz;
    channel_labels = opts.channel_labels;

    % If the length of channel_labels does not match the number of channels
    % in the recording, relabel channels as their index
    if length(channel_labels) ~= size(data_tbl.segment_data{1},1)
        channel_labels = num2cell(1:size(data_tbl.segment_data{1},1));
    end

    if ~exist(sprintf("%s/%s/", save_fig_loc,imprint_changes_label), 'dir')
        mkdir(sprintf("%s/%s/", save_fig_loc,imprint_changes_label))
    end

    if size(cell_imprint,1) > 1 & ~isnan(sz)
        fprintf("Can only use opts.sz if only one seizure is loaded into function \n")
        return
    end

    for s = 1:size(cell_imprint,1)
        f = figure('Visible', fig_visible);
        f.Position = [100,100, 1600, 1000];

        fs = data_tbl.segment_fs(s);
%         subplot(1,9,1:8)
        eeg_dat = data_tbl.segment_data{s,1};
        pre_dat = eeg_dat(:,1:109*fs);

        offset = round((max(max(pre_dat))-min(min(pre_dat)))/3);
        % Extract removed segments (this is where abnormal activity was
        % detected)
        preict_mad = cell_imprint.cell_pre_features_mad{s,1};
        abnorm_act = isnan(preict_mad(:, 1:(end-1)));
        abnorm_act_eeg = repelem(abnorm_act, 1, 512/(8));
        abnorm_act_bin = repelem(abnorm_act, 1, 512/(8^2));

        % Extract segments of EEG which have abnormal preictal activity
        % detected
        abnorm_eeg = pre_dat.*abnorm_act_eeg;
        abnorm_eeg(abnorm_eeg == 0) = NaN;

        hold on
        vis_eeg(pre_dat,...%(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", channel_labels,...
            "PlotNewFig", false, "Color", [0.6,0.6,0.6], "Offset", offset);
        vis_eeg(abnorm_eeg,... %(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", channel_labels,...
            "PlotNewFig", false, "Color", [0.7,0,0], "Offset", offset);
        hold off
        colormap([1,1,1;1,0.2,0])
        set(gca, "XTick", 0:80:880, "XTickLabel", -120:10:-10)
        
        if isnan(sz)
            sz = s;
            sz_type = data_tbl.ilae_sz_type(sz);
        else
            sz_type = data_tbl.ilae_sz_type;
        end
        
        sgtitle(sprintf("%s seizure %d (%s)", string(data_tbl(s,:).patient_id), sz, sz_type))
        if save_fig == 1
            saveas(f,sprintf("%s/%s/preict_sz%d.png", save_fig_loc, ...
                imprint_changes_label, sz))
        end
        
        
    end