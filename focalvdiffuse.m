data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];
patient = "1005";
% Load json data to get xyz coordianes of onset channels
pat_onset = final_output(final_output.Patient_id == patient,:);
filelist = dir(fullfile(strcat(path_pipeline_exports, "/UCLH", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  % remove folders from list
folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder

json_channels = json_data(1).channel_details;
json_chan_name = strrep(json_channels.chan_name, ' ','');

incl_chan = ismember(json_chan_name, strrep(pat_onset.channel_names{1,1},' ',''));
included_channels = json_channels(incl_chan,:);
incl_chan_loc = cat(2,included_channels.location_orig{:})';

onset_xyz = incl_chan_loc(included_channels.is_soz,:);

% Compute Euclidean distances between electrodes
dists = squareform(pdist(onset_xyz));
% Extract the maximum distance
max(max(dists))

% Need to determime a threshold to suggest if onset is focal or diffuse
% I think this is on the same scale as the distances in the atlas so I
% could divide by the maximum there and see if there appears to be a
% separating line

%% Decided to do this on a region level instead
patients = final_output.Patient_id;
for pat = 1:length(patients)
    pat_onset = final_output(final_output.Patient_id == patients(pat),:);
    pat_roi = string(pat_onset.ROI_ids{:});
    onset_roi =  pat_roi(find(pat_onset.labelled_onset_roi{:}));
    pat_dists = dists125(find(ismember(string(atlas.name{3}), onset_roi)),...
        find(ismember(string(atlas.name{3}), onset_roi)));
    final_output.max_clo_dist(pat) = max(max(pat_dists));
    sz_max_dist = nan(length(pat_onset.Segment_ids{:}),1);
    for sz = 1:length(pat_onset.Segment_ids{:})
        if sum(pat_onset.imprint_roi{:}(:,sz)) ~= 0
            sz_roi =  pat_roi(find(pat_onset.imprint_roi{:}(:,sz)));
            sz_dists = dists125(find(ismember(string(atlas.name{3}), sz_roi)),...
                find(ismember(string(atlas.name{3}), sz_roi)));
            sz_max_dist(sz) = max(max(sz_dists));
        else
            continue
        end
    end
    final_output.max_imprint_dists(pat) = {sz_max_dist};
end