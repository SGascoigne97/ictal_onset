load('roi_info/ATLAS.mat')
load('roi_info/lhSurf.mat') 
load('roi_info/rhSurf.mat')
xyz = atlas.xyz{3,1};

%% Load patient data (pre-processed)
data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal';
path_pipeline_exports = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/export_ictal';
patient = 'UCLH851';
load(sprintf('%s/%s.mat', data_location, patient))
%% load json_data
filelist = dir(fullfile(strcat(path_pipeline_exports, "/", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  % remove folders from list
folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder

%% Extract resected region and check that it matches with PT
channel_details = json_data(1).channel_details;
recorded_channels = channel_details(ismember(channel_details.chan_name, data_export.segment_channel_labels{1}),:);

recorded_roi_all_atlas = cat(2,recorded_channels.ROIname{:});
recorded_regions = unique(recorded_roi_all_atlas(3,:)');

resected_regions = unique(recorded_roi_all_atlas(3,find(recorded_channels.is_resected5))');

recorded_xyz = xyz(find(ismember(atlas.name{3,1}, recorded_regions)),:);
resected_xyz = xyz(find(ismember(atlas.name{3,1}, resected_regions)),:);

channelx = [];
channely = [];
channelz = [];
res_channelx = [];
res_channely = [];
res_channelz = [];

for i = 1:height(recorded_channels)
    location = recorded_channels.location_orig{i};
    if ~recorded_channels.is_resected3(i)
        channelx(end+1) = location(1);
        channely(end+1) = location(2);
        channelz(end+1) = location(3);
    else
        res_channelx(end+1) = location(1);
        res_channely(end+1) = location(2);
        res_channelz(end+1) = location(3);
    end

end


% Plot channel and region locations
figure()
scatter3(recorded_xyz(:,1),recorded_xyz(:,2),recorded_xyz(:,3), 'black')
hold on
Hl = patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
% camlight('headlight','infinite')
Hv = patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
% camlight('headlight','infinite')
scatter3(resected_xyz(:,1),resected_xyz(:,2), resected_xyz(:,3),'red', 'filled')
scatter3(res_channelx, res_channely, res_channelz, 'green', 'filled');
axis equal
hold off
