% Comparing detected onset with clinically labelled onset 
%load("tables/all_pat_with_focal_2.mat");

%%
sz_count_thresh = 0.01; 

    in_min_sz_tab = table();

     cm = [0.9,0.9,0.9; % neither (0:grey)
        0,0,1; % CLO not captured (1: blue)
        1,1,0; % in imprint, not clo (2:yellow)
         0,1,0];% CLO captured (3:yellow)
    cm = interp1(cm, 1:0.01:size(cm,1));
%
for pat = 1:size(final_output,1)
    patient = final_output.Patient_id{pat};
    pat_onset = final_output(string(final_output.Patient_id) == patient,: );
    if all(isnan(pat_onset.clo_chan{:}))
        fprintf("Patient %s does not have clinically labelled onset \n", patient)
        continue
    end
    %%
%     figure(pat)
%     subplot(2,4,1:3)
%     imagesc(2*pat_onset.imprint_chan{:}+pat_onset.clo_chan{:})
%     clim([0 3])
%     set(gca, 'YTick', 1:length(pat_onset.channel_names{:}), 'YTickLabels', pat_onset.channel_names{:} )
%     title("Channel-wise")
%     subplot(2,4,4)
%     imagesc(2*pat_onset.imprint_chan{:}+pat_onset.clo_chan{:})
%     clim([0 3])
%     set(gca, 'YTick', 1:length(pat_onset.channel_names{:}), 'YTickLabels', pat_onset.channel_names{:} )
% % 
%     subplot(2,4,5:7)
%     imagesc(2*pat_onset.imprint_roi{:}+pat_onset.clo_roi{:})
%     clim([0 3])
%     set(gca, 'YTick', 1:length(pat_onset.roi_names{:}), 'YTickLabel', pat_onset.roi_names{:} )
%     title("Region-wise")
%     sgtitle(sprintf("%s comparing onsets (CLO and imprint)", patient))


    % In one sz comparison
    min_sz_count = ceil(sz_count_thresh*length(pat_onset.Segment_ids{:}));
    in_min_sz_chan = sum(pat_onset.imprint_chan{:},2)>=min_sz_count;
    in_min_sz_roi = sum(pat_onset.imprint_roi{:},2)>=min_sz_count;

%     subplot(2,4,4)
%     imagesc(2*in_min_sz_chan+pat_onset.clo_chan{:})
%     clim([0 3])
%     set(gca, 'YTick', [], 'XTick',[])
%     subplot(2,4,8)
%     imagesc(2*in_min_sz_roi+pat_onset.clo_roi{:})
%     clim([0 3])
%     set(gca, 'YTick', [], 'XTick',[])
%     roi_names = pat_onset.roi_names{:};
%     roi_names = strrep(roi_names, 'r.', 'ctx-rh-');
%     roi_names = strrep(roi_names, 'l.', 'ctx-lh-');

    %plotBrain(roi_names, 2*in_min_sz_roi+pat_onset.clo_roi{:},...
      %  cm, 'atlas', 'lausanne120_aseg')

% 
    %
    pat_comp_tab = table();
    for sz = 1:size(pat_onset.imprint_chan{:},2)
        pat_comp_tab.patient_id(sz) = {patient};
        pat_comp_tab.jac_chan(sz) = jaccard(pat_onset.clo_chan{:},logical(pat_onset.imprint_chan{:}(:,sz)));
        pat_comp_tab.jac_roi(sz) = jaccard(logical(pat_onset.clo_roi{:}),logical(pat_onset.imprint_roi{:}(:,sz)));
        pat_comp_tab.perc_clo_chan(sz) = sum(pat_onset.imprint_chan{:}(:,sz)+pat_onset.clo_chan{:} == 2)/sum(pat_onset.clo_chan{:});
        pat_comp_tab.perc_clo_roi(sz) =  sum(pat_onset.imprint_roi{:}(:,sz)+pat_onset.clo_roi{:} == 2)/sum(pat_onset.clo_roi{:});
    end
    
    if exist('all_pat_comp', 'var')
        all_pat_comp = [all_pat_comp; pat_comp_tab];
    else
        all_pat_comp = pat_comp_tab;
    end

    pat_surg_out = pat_onset.Surgery_outcome{:};
    pat_surg_year = pat_onset.("Surgery year");
    pat_out_year = pat_onset.("Outcome year"){:};
    year_1 = pat_out_year - pat_surg_year == 1;
    

    for sz = 1:size(pat_onset.imprint_chan{:},2)
        in_min_sz_tab.patient_id(pat) = {patient};
        in_min_sz_tab.jac_chan(pat) = jaccard(pat_onset.clo_chan{:},in_min_sz_chan);
        in_min_sz_tab.jac_roi(pat) = jaccard(logical(pat_onset.clo_roi{:}),in_min_sz_roi);
        in_min_sz_tab.perc_clo_chan(pat) = sum(in_min_sz_chan+pat_onset.clo_chan{:} == 2)/sum(pat_onset.clo_chan{:});
        in_min_sz_tab.perc_clo_roi(pat) =  sum(in_min_sz_roi+pat_onset.clo_roi{:} == 2)/sum(pat_onset.clo_roi{:});
        in_min_sz_tab.outcome(pat) = pat_surg_out(year_1);
    end
    

end

%
in_min_sz_tab = in_min_sz_tab(~cellfun(@isempty,in_min_sz_tab.patient_id),:);

figure()
subplot(251)
boxchart(ones(1,size(in_min_sz_tab,1)), in_min_sz_tab.perc_clo_chan)
ylabel("Percentage of CLO captured")
hold on
swarmchart(ones(1,size(in_min_sz_tab,1)), in_min_sz_tab.perc_clo_chan,'filled')
hold off
title("Percentage clo channels detected")
subplot(2,5,2:5)
boxchart(in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).outcome, in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).perc_clo_chan)
xlabel("Surgical outcome")
hold on
swarmchart(in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).outcome, in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).perc_clo_chan,'filled')
hold off

subplot(256)
boxchart(ones(1,size(in_min_sz_tab,1)), in_min_sz_tab.perc_clo_roi)
ylabel("Percentage of CLO captured")
hold on
swarmchart(ones(1,size(in_min_sz_tab,1)), in_min_sz_tab.perc_clo_roi,'filled')
hold off
title("Percentage clo regions detected")
subplot(2,5,7:10)
boxchart(in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).outcome, in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).perc_clo_roi)
xlabel("Surgical outcome")
hold on
swarmchart(in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).outcome, in_min_sz_tab(in_min_sz_tab.outcome ~= 8,:).perc_clo_roi,'filled')
hold off
sgtitle("Percentage of clo in min sz (>25%)")
