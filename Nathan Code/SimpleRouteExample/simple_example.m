% hypothesis: Seizures should follow white matter, when they stop doing
% this they have spread outside the implanted area

addpath(genpath('Brewermap'))
%% load atlas
load('/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/roi_info/ATLAS.mat', 'atlas');
load('/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/roi_info/lhSurf.mat', 'lf', 'lv');
load('/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/roi_info/rhSurf.mat', 'rf', 'rv');

%% load imprint code
%imprint_location = '/Users/b9056471/Documents/GitHub/ictal_onset/';
%addpath(genpath(imprint_location));

%% load connectivity matricies
load('sFrances.mat', 's');
s_table = struct2table(s);

%% load touching matrix 36 60 125 250
%% scale 125 touching matrix left hand side not symetric???????
scale = 60;
[names, touching, roiIndex] = load_touching_matrix(scale);

%% everybodies favourite patient
patient = 'UCLH1005';
conData = s_table.scale60(strcmp(s_table.subjID,patient(5:end)));

%% load patient data
load("1005_data_export.mat", 'data_export');
load("1005_json_data.mat", 'json_data');
   
%% calcaulate imprints
[pat_data, ~, cell_imprint,  sz_count_pat] = calc_imprint(data_export, data_export, 'window_size', 1, 'min_sz_count', 1, 'window_overlap', 7/8,'folder','onset_calcs_8');
%% convert to regions
cell_imprint_regions = imprint_to_regions(cell_imprint, json_data, data_export);

%% create white matter touching matrix
wm_names = conData{1}.names;
wm_names = regexp(wm_names,"[^']*",'match','once');
wm_touching = logical(conData{1}.adjCount);

%% match touching matrix to wm touching matrix
idx = find_preserved(names,wm_names);
idx2 = find_preserved(wm_names,names);
touching_wm_size = zeros(size(wm_touching));
touching_wm_size(idx2,idx2) = touching(idx,idx);

region_names = cell_imprint_regions{1,"region_names_" + string(roiIndex)};
region_indexes = find_preserved(wm_names,region_names);

%% loop over every seizure
for seg = 1:height(cell_imprint_regions)
    region_imprint = logical(cell_imprint_regions{seg,"mean_imprint_region_" + string(roiIndex)}{1});
    t0 = find((any(region_imprint,1)),1);
    if size(region_imprint,2) <= t0 + 8
        onset_region_imprint = region_imprint(:,t0:end);
    else
        onset_region_imprint = region_imprint(:,t0:t0+8);
    end
    onset_region_imprint_no_duplicates = onset_region_imprint(any(onset_region_imprint,1),:);
    onset_region_imprint_no_duplicates = onset_region_imprint_no_duplicates(:,[true, any(onset_region_imprint_no_duplicates(:,1:end-1) ~= onset_region_imprint_no_duplicates(:,2:end))]);
    if size(onset_region_imprint_no_duplicates,2) > 1
        seizureRoute_onset = showRoute(wm_names, onset_region_imprint, touching_wm_size, wm_touching, region_indexes, roiIndex, atlas, lf, lv, rf, rv);
        seizureRoute_all = showRoute(wm_names, region_imprint, touching_wm_size, wm_touching, region_indexes, roiIndex, atlas, lf, lv, rf, rv);
    end
end
