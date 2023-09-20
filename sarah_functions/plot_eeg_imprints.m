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
        f = figure('Visible', 'on');
        f.Position = [100,100, 1600, 1000];
%         subplot(1,9,1:8)
        eeg_dat = data_tbl.segment_data{sz,1};

        offset = round((max(max(eeg_dat))-min(min(eeg_dat)))/3);
        imprint = cell_imprint.cell_imprint{sz,1};
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
   
        onset_heatmap = zeros(size(imprint_bin,1), size(imprint_bin,2));
        onset_heatmap(:,find(sum(imprint_bin,1)>0, 1, 'first') +[0:7]) = ...
            imprint_bin(:,find(sum(imprint_bin,1)>0, 1, 'first') +[0:7]);

        hold on
        imagesc(imprint_bin, 'AlphaData', 0.5)
        imagesc(2*onset_heatmap, "AlphaData", 0.7)
        vis_eeg(eeg_dat,...%(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", data_tbl.segment_channel_labels{1,1}(:),...
            "PlotNewFig", false, "Color", [0.2,0.2,0.2], "Offset", offset);
        vis_eeg(imprint_eeg_dat,... %(:,(512*60):(512*(180+data_tbl.duration(sz)))),...
            512/8, "ChannelNames", data_tbl.segment_channel_labels{1,1},...
            "PlotNewFig", false, "Color", [0,0,0.7], "Offset", offset);
        hold off
        colormap([1,1,1;0,0.2,1;0,0.5,0])
        clim([0,2])
        xline(120*8, LineWidth=2, Color="red")
        xline(find(sum(imprint_bin),1, 'first'), LineWidth=2, Color="blue")
        xline((120+data_tbl.duration(sz))*8, LineWidth=2, Color="red")
        set(gca, "XTick", 0:10*8:size(eeg_dat,2)/8/8, "XTickLabel", -120+(0:10:size(eeg_dat,2)/8^3))
        xlim([90*8, max([150+find(sum(imprint_bin),1, 'first')/(512/8), 120+data_tbl.duration(sz)])*8])
        
%         subplot(1,9,9)
%         imagesc(flipud(2*double(sum(imprint_bin(:,find(sum(imprint_bin),1, 'first')+[0:7]),2)>0)), 'AlphaData', 0.5)
%         axis off
%         clim([0,2])
        
        if isnan(opts.sz)
            sgtitle(sprintf("%s seizure %d (%s)", string(data_tbl(sz,:).patient_id), sz, data_tbl.ilae_sz_type(sz)))
            saveas(f,sprintf("%s/%s/sz%d.png", opts.save_fig_loc, ...
                opts.imprint_changes_label, sz))
        else
            sz_type = data_tbl.ilae_sz_type;
            if ismissing(sz_type)
                sz_type = "NaN";
            end
            sgtitle(sprintf("%s seizure %d (%s)", string(data_tbl.patient_id), opts.sz, sz_type))
            saveas(f,sprintf("%s/%s/sz%d.png", opts.save_fig_loc, ...
                opts.imprint_changes_label, opts.sz))
        end
        
    end