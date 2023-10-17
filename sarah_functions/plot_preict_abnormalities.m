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
    end

    if ~exist(sprintf("%s/%s/", opts.save_fig_loc,opts.imprint_changes_label), 'dir')
        mkdir(sprintf("%s/%s/", opts.save_fig_loc,opts.imprint_changes_label))
    end

    if size(cell_imprint,1) > 1 & ~isnan(opts.sz)
        fprintf("Can only use opts.sz if only one seizure is loaded into function \n")
        return
    end

    for sz = 1:size(cell_imprint,1)
        f = figure('Visible', opts.fig_visible);
        f.Position = [100,100, 1600, 1000];

        fs = data_tbl.segment_fs(sz);
%         subplot(1,9,1:8)
        eeg_dat = data_tbl.segment_data{sz,1};
        pre_dat = eeg_dat(:,1:109*fs);

        offset = round((max(max(pre_dat))-min(min(pre_dat)))/3);
        % Extract removed segments (this is where abnormal activity was
        % detected)
        preict_mad = cell_imprint.cell_pre_features_mad{sz,1};
        abnorm_act = isnan(preict_mad(:, 1:(end-1)));
        abnorm_act_eeg = repelem(abnorm_act, 1, 512/(8));
        abnorm_act_bin = repelem(abnorm_act, 1, 512/(8^2));

        % Extract segments of EEG which have abnormal preictal activity
        % detected
        abnorm_eeg = pre_dat.*abnorm_act_eeg;
        abnorm_eeg(abnorm_eeg == 0) = NaN;

        hold on
        vis_eeg(pre_dat,...%(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", data_tbl.segment_channel_labels{1,1}(:),...
            "PlotNewFig", false, "Color", [0.6,0.6,0.6], "Offset", offset);
        vis_eeg(abnorm_eeg,... %(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", data_tbl.segment_channel_labels{1,1},...
            "PlotNewFig", false, "Color", [0.7,0,0], "Offset", offset);
        hold off
        colormap([1,1,1;1,0.2,0])
        set(gca, "XTick", 0:80:880, "XTickLabel", -120:10:-10)
        
        if isnan(opts.sz)
            sgtitle(sprintf("%s seizure %d (%s)", string(data_tbl(sz,:).patient_id), sz, data_tbl.ilae_sz_type(sz)))
            saveas(f,sprintf("%s/%s/preict_sz%d.png", opts.save_fig_loc, ...
                opts.imprint_changes_label, sz))
        else
            sz_type = data_tbl.ilae_sz_type;
            if ismissing(sz_type)
                sz_type = "NaN";
            end
            sgtitle(sprintf("%s seizure %d (%s)", string(data_tbl.patient_id), opts.sz, sz_type))
            saveas(f,sprintf("%s/%s/preict_sz%d.png", opts.save_fig_loc, ...
                opts.imprint_changes_label, opts.sz))
        end
        
    end