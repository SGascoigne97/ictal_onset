
load('roi_info/ATLAS_.mat')
load('roi_info/retains.mat')
load('roi_info/atlasinfo.mat')
load('roi_info/ATLAS.mat')
%load('roi_info/lhSurf.mat') 
%load('roi_info/rhSurf.mat')

addpath('sarah_functions/')
addpath(genpath('help_functions/'))
addpath('Simple-Brain-Plot-main/')


%%
load("final_output.mat");
    
for pat = 1:length(final_output.Patient_id)
    %%
    patient = final_output.Patient_id(pat);
    pat_onset = final_output(final_output.Patient_id == patient,:);
    
    % xyz coordinates for regions in 125 atlas
    % 
    %xyz = atlas.xyz{3,1};
    roi_names = atlas.name{3,1};
    
    patient_rois = pat_onset.ROI_ids{:};
    incl_roi = ismember(roi_names, patient_rois);
    %patient_xyz = xyz(incl_roi,:); % Careful - the order of regions is not preserved here 
    
    labelled_onset = pat_onset.labelled_onset_roi{1,1};
    labelled_onset(isnan(labelled_onset)) = 0;
    onset_names = patient_rois(find(labelled_onset));
    resected = pat_onset.resected_roi{1,1};
    resected(isnan(resected)) = 0;
    resected_names = patient_rois(find(resected));
    onset_and_resec = labelled_onset + resected == 2;
    o_and_r_names = patient_rois(find(onset_and_resec));
     
    % We want to look at the json data to see which regions were resected including all channels (before preprocessing)
    
    data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
    path_pipeline_exports = [data_location, 'export_ictal'];
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/UCLH", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
    
    
    channels_before_preproc = json_data(1).channel_details;
    channels_before_preproc = channels_before_preproc(~cellfun(...
        'isempty', channels_before_preproc.ROIids),:);
    roi_names = cat(2, channels_before_preproc.ROIname{:});
    roi_60 = roi_names(2,:);
    unq_roi = unique(roi_60);
    
    resected_roi = chan_to_roi_crit(channels_before_preproc.is_resected5, roi_60, "threshold", 0.25);
    
    channels_retained = ismember(cat(2,channels_before_preproc.chan_name),...
        pat_onset.channel_names{:});
    onset_roi = chan_to_roi_crit(pat_onset.labelled_onset_chan{:}, ...
        roi_60(find(channels_retained)), "threshold", 0.25);

    resected_roi(isnan(resected_roi)) = 0;
    onset_roi(isnan(onset_roi)) = 0;
    
    unq_roi = strrep(unq_roi, 'r.', 'ctx-rh-');
    unq_roi = strrep(unq_roi, 'l.', 'ctx-lh-');
    resected_names = unq_roi(find(resected_roi));
    pat_roi = unique(roi_names(2,find(channels_retained)));
    pat_roi = strrep(pat_roi, 'r.', 'ctx-rh-');
    pat_roi = strrep(pat_roi, 'l.', 'ctx-lh-');
    onset_names = pat_roi(find(onset_roi));
    
    
    cm = [0.9,0.9,0.9; %  recorded but no channels remaining (-1: pale grey)
        0.5,0.5,0.5; % neither (0:grey)
        1,1,0; % +onset- resec (1:yellow)
        0,0,1; % -onset +resec (2:blue)
        0,1,0];% +onset +resec (3:green)
    cm = interp1(cm, 1:0.01:size(cm,1));

    not_recorded = ~ismember(unq_roi, pat_roi);
    onset = ismember(unq_roi, onset_names);
    col_sch = onset + 2*resected_roi ;
    col_sch(find(not_recorded)) = -1;
    
    plotBrain(unq_roi, col_sch,cm, 'atlas', 'lausanne120_aseg', 'limits', [-1,3], ...
        'savePath', char(sprintf('./figures/onset_resec/%s', patient)))
end

