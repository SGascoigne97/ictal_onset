function plot_eeg_imprints(data_tbl, cell_imprint, opts)
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
        opts.imprint_changes_label = "original" % Use this if you are 
        % making changes to the algorithm and would like to track changes to imprints and onsets 
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

    if ~exist(sprintf("%s/%s/", save_fig_loc, imprint_changes_label), 'dir')
        mkdir(sprintf("%s/%s/", save_fig_loc, imprint_changes_label))
    end

    if size(cell_imprint,1) > 1 & ~isnan(sz)
        fprintf("Can only use opts.sz if only one seizure is loaded into function \n")
        return
    end

    for s = 1:size(cell_imprint,1)
        f = figure('Visible', fig_visible);
        f.Position = [100,100, 1600, 1000];
%         subplot(1,9,1:8)
        eeg_dat = data_tbl.segment_data{s,1};

        offset = round((max(max(eeg_dat))-min(min(eeg_dat)))/3);
        imprint = cell_imprint.cell_imprint{s,1};
        imprint_bin = [zeros(size(imprint,1),110*8), imprint];
        imprint_bin = repelem(imprint_bin,offset,1);
        for chan = 1:size(eeg_dat,1)
            top = 900 + (chan-1)*offset;
            bottom = chan*offset;
            imprint_bin(top:bottom,:) = 0;
        end
        imprint_bin = [imprint_bin; zeros(round(offset/2),size(imprint_bin,2))];
        imprint_bin = imprint_bin(round(offset/2):end,:);
    
        % Extract segments of EEG which are included in imprint
        % Add preictal and ictal segments to imprint
        imprint_full = [zeros(size(imprint,1),110*8), imprint, zeros(size(imprint,1),120*8)];
        imprint_eeg = repelem(imprint_full, 1, round(size(eeg_dat,2)/size(imprint_full,2)));
        % Add in zeros at each side so the sizes are equivalent 
        dif = size(eeg_dat,2) - size(imprint_eeg,2);
        if dif > 0
            imprint_eeg = [zeros(size(imprint_eeg,1), round(dif/2)), imprint_eeg, zeros(size(imprint_eeg,1), round(dif/2)-1)];
        elseif dif < 0
            dif = abs(dif);
            imprint_eeg = imprint_eeg(:, round(dif/2):end-round(dif/2));
        end
    
        imprint_eeg(imprint_eeg == 0) = NaN;
        imprint_eeg_dat = eeg_dat.*imprint_eeg;
   
        if ~isempty(find(sum(imprint_bin,1)>0, 1, 'first'))
            onset_heatmap = zeros(size(imprint_bin,1), size(imprint_bin,2));
            onset_heatmap(:,find(sum(imprint_bin,1)>0, 1, 'first') +[0:7]) = ...
                imprint_bin(:,find(sum(imprint_bin,1)>0, 1, 'first') +[0:7]);
        else 
            continue
        end

        hold on
        imagesc(imprint_bin, 'AlphaData', 0.5)
        imagesc(2*onset_heatmap, "AlphaData", 0.7)
        vis_eeg(eeg_dat,...%(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", channel_labels,...
            "PlotNewFig", false, "Color", [0.2,0.2,0.2], "Offset", offset);
        vis_eeg(imprint_eeg_dat,... %(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", channel_labels,...
            "PlotNewFig", false, "Color", [0,0,0.7], "Offset", offset);
        hold off
        colormap([1,1,1;0,0.2,1;0,0.5,0])
        clim([0,2])
        xline(120*8, LineWidth=2, Color="red")
        xline(find(sum(imprint_bin),1, 'first'), LineWidth=2, Color="blue")
        xline((120+data_tbl.duration(s))*8, LineWidth=2, Color="red")
        set(gca, "XTick", 0:10*8:size(eeg_dat,2)/8/8, "XTickLabel", -120+(0:10:size(eeg_dat,2)/8^3))
        xlim([90*8, max([150+find(sum(imprint_bin),1, 'first')/(512/8), 120+data_tbl.duration(s)])*8])


        if isnan(sz)
            sz = s;
            sz_type = data_tbl.ilae_sz_type(sz);
        else
            sz_type = data_tbl.ilae_sz_type;
        end
        
        sgtitle(sprintf("%s seizure %d (%s)", string(data_tbl(s,:).patient_id), sz, sz_type))
        if save_fig == 1
            saveas(f,sprintf("%s/%s/sz%d.png", save_fig_loc, ...
                    imprint_changes_label, sz))
        end
       
        
    end