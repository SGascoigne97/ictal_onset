pat =  16; 
pat_onset = final_output(pat,:);

unq_roi = pat_onset.roi_names_120{:};
unq_roi = strrep(unq_roi, 'r.', 'ctx-rh-');
unq_roi = strrep(unq_roi, 'l.', 'ctx-lh-');

onset = pat_onset.imprint_roi_120{:};

onset_across = double((sum(onset,2)/size(onset,2)) >= 0.5);
cm = [0.7,0.7,0.7; 0,0,1];
plotBrain(unq_roi, onset_across, cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1],...
    'savePath', char(sprintf('./figures/onset_var/%s_imprint_across', pat_onset.Patient_id{:})))
cm = [0.7,0.7,0.7; 1,0,0];
plotBrain(unq_roi, pat_onset.clo_roi_120{:}, cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1],...
    'savePath', char(sprintf('./figures/onset_var/%s_clo', pat_onset.Patient_id{:})))
cm = [0.7,0.7,0.7; 0,1,0];
plotBrain(unq_roi, pat_onset.resected_roi_120{:}, cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1],...
    'savePath', char(sprintf('./figures/onset_var/%s_resec', pat_onset.Patient_id{:})))

% Plot proportion of seizures in onset for each region
cm = [0.7,0.7,0.7; 0,0,1; 0,0,1];
cm = interp1(cm, 1:0.01:size(cm,1));
prop_sz = sum(onset,2)/size(onset,2);
plotBrain(unq_roi, prop_sz, cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1],...
    'savePath', char(sprintf('./figures/onset_var/%s_prop_sz', pat_onset.Patient_id{:})))

%, ...
   % 'savePath', char(sprintf('./figures/onset_var/%s', patient)))