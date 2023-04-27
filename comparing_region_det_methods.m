data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];
load('roi_info/lhSurf.mat') 
load('roi_info/rhSurf.mat')

%% Plot resections on region and channel level for all patients
for pat = 1:length(patients)
    patient = patients(pat);
    pat_onset = final_output(final_output.Patient_id == patient,:);
    pat_resec_roi = pat_onset.resected_roi{:};
    pat_resec_roi(isnan(pat_resec_roi)) = 0;
    pat_resec_chan = pat_onset.resected_chan{:};
    
    pat_roi_names = pat_onset.ROI_ids{:};
    pat_resec_roi_names = pat_roi_names(find(pat_resec_roi));
    
    roi_atlas = ismember(atlas.name{3,1}, pat_roi_names);
    roi_resec_atlas = ismember(atlas.name{3,1}, pat_resec_roi_names);
    
    xyz = atlas.xyz{3,1};
    
    roi_xyz = xyz(roi_atlas,:);
    roi_resec_xyz = xyz(roi_resec_atlas,:);
    
    %% Load json data to obtain channel xyz location
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/UCLH", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder

    % Remove channels with no ROI mapping    
    empty_roi = cellfun(@isempty, json_data(1).channel_details.ROIname);
    chan_with_roi = json_data(1).channel_details(~empty_roi,:);

    pat_channel_details = chan_with_roi(ismember(...
        strrep(string(chan_with_roi.chan_name), ' ', ''),...
        strrep(string(pat_onset.channel_names{:}), ' ', '')),:);
    pat_chan_locations = cat(2, pat_channel_details.location_orig{:})';

    pat_chan_resec = pat_onset.resected_chan{:};
    pat_chan_resec = pat_chan_resec(ismember(string(pat_channel_details.chan_name),...
        string(pat_onset.channel_names{:})));

    pat_resec_chan_xyz = pat_chan_locations(find(pat_chan_resec),:);
    
    %% Plot on the same surface
    figure(pat)
    scatter3(roi_xyz(:,1),roi_xyz(:,2),roi_xyz(:,3),'yellow')
    hold on
    Hl = patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    Hv = patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    scatter3(roi_resec_xyz(:,1),roi_resec_xyz(:,2),roi_resec_xyz(:,3),'red', 'filled')
    scatter3(pat_resec_chan_xyz(:,1),pat_resec_chan_xyz(:,2),pat_resec_chan_xyz(:,3),'blue', 'filled')
    hold off
    title(sprintf("Patient %s resected regions (red) and channels (blue)", patient))
    saveas(gcf,sprintf('figures/check_locations/%s', patient), 'png')
end

 

    %% Plot resections on region and channel level for all patients
for pat = 1:length(patients)
    patient = patients(pat);
    pat_onset = final_output(final_output.Patient_id == patient,:);
    pat_resec_roi = pat_onset.resected_roi{:};
    pat_resec_roi(isnan(pat_resec_roi)) = 0;
    pat_resec_chan = pat_onset.resected_chan{:};
    
    pat_roi_names = pat_onset.ROI_ids{:};
    pat_resec_roi_names = pat_roi_names(find(pat_resec_roi));
    
    roi_atlas = ismember(atlas.name{3,1}, pat_roi_names);
    roi_resec_atlas = ismember(atlas.name{3,1}, pat_resec_roi_names);
    
    xyz = atlas.xyz{3,1};
    
    roi_xyz = xyz(roi_atlas,:);
    roi_resec_xyz = xyz(roi_resec_atlas,:);

    % Compare old and new methods of labelling resected regions
    chan_names = pat_onset.channel_names{:};
    resec_chan = chan_names(find(pat_onset.resected_chan{:}));
    unq_roi = pat_onset.ROI_ids;
    %% Load json data to obtain channel xyz location
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/UCLH", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder

    empty_roi = cellfun(@isempty, json_data(1).channel_details.ROIname);
    chan_with_roi = json_data(1).channel_details(~empty_roi,:);

    pat_channel_details = chan_with_roi(ismember(...
        strrep(string(chan_with_roi.chan_name), ' ', ''),...
        strrep(string(pat_onset.channel_names{:}), ' ', '')),:);

    % Remove channels that are missing ROI information
    chan_roi_names_all = cat(2, pat_channel_details.ROIname{:});
    chan_roi_names = string(chan_roi_names_all(3,:))';

    unq_roi = unique(chan_roi_names);

    pat_chan_resec = pat_onset.resected_chan{:};
    pat_chan_resec = pat_chan_resec(ismember(string(pat_channel_details.chan_name),...
        string(pat_onset.channel_names{:})));
    
    roi_resec_names = chan_roi_names(find(pat_chan_resec));

    chan_in_roi = zeros(1, length(unq_roi));
    for grp = 1:length(unq_roi)
        chan_in_roi(grp) = sum(ismember(unq_roi(grp), roi_resec_names)) >0;
    end

    roi_atlas = ismember(atlas.name{3,1}, pat_roi_names);
    old_roi_resec_atlas = ismember(atlas.name{3,1}, unq_roi(find(chan_in_roi)));
    
    xyz = atlas.xyz{3,1};
    old_roi_resec_xyz = xyz(old_roi_resec_atlas,:);


    %% Plot on the same surface
    figure(pat)
    scatter3(roi_xyz(:,1),roi_xyz(:,2),roi_xyz(:,3),'yellow')
    hold on
    Hl = patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    Hv = patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
    scatter3(roi_resec_xyz(:,1),roi_resec_xyz(:,2),roi_resec_xyz(:,3),'red', 'filled')
    scatter3(old_roi_resec_xyz(:,1)+0.5,old_roi_resec_xyz(:,2)+0.5,old_roi_resec_xyz(:,3)+0.5,'blue', 'filled')
    title(sprintf("Patient %s resected regions using old (blue) and new (red) code", patient))
    saveas(gcf,sprintf('figures/check_methods/%s', patient), 'png')
    hold off
end