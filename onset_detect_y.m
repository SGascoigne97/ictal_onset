clear all
close all
%%
data_location = '~/Desktop/Sarah_onset/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];

% Add paths to all functions required
addpath(genpath('help_functions'))
addpath('lib_biomarkers')
addpath('lib_dataflow')
addpath('sarah_functions')
addpath(genpath('../ieeg-norm-map-pipeline/lib/'))

%%
% List patients with pre-processed data
patients_dir = dir(path_pipeline_exports);
patients = {patients_dir(4:end).name};

% Choose atlas
atl = 3;%1=scale36,2=60,3=125,4=250
% Choose threshold for allowing propagated activity (seconds)
det = 0;

min_sz = 8;
chan_to_roi_thresh = 0.25; % One channel in region is sufficient to include region

onset_calc_loc = "onset_calcs"; % Specify folder to store imprint values in
% Need a new folder if using a different subset of the data/different data
% as it will load previous save if folder is not empty

% For each patient, compute onset based on imprint
pat = 1%for patients
patient = patients{pat};

    %if exist(sprintf('%s/%s.mat', data_location, patient), 'file')
        load(sprintf('%s/%s.mat', data_location, patient));
%     %else
%         fprintf('Patient %s does not have saved data \n', patient)
%         continue 
%     end
    pat_data = data_export;
    
    % load json_data
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
    
    % Apply inclusion criteria
    [pat_data, ~, ~] = incl_crit(pat_data, 'min_sz_count', min_sz);
    % Append patient data to data table (if inclusion criteria are met)
%     if size(pat_data,1) < min_sz
%         fprintf('Patient %s does not meet inclusion criteria \n', patient)
%         continue
%     end
    pat_meta = pat_data(:,1:(end-1));

    % Remove seizures from json data that do not meet inclusion criteria
    json_data = json_data(contains(string(extractfield(cat(2,json_data.x_id),...
        'x_oid')), string(pat_data.segment_id)));
    
    % Compute onsets
    mkdir(onset_calc_loc)

%     if size(pat_data,1) > 0
        % Compute imprint for all seizures
        [pat_data, pat_meta, cell_imprint,  sz_count_pat] = calc_imprint(pat_data, pat_meta);
%     else 
%         fprintf("Patient %s does not meet inclusion criteria \n", patient)
%         continue
%     end
%     
%     % Check order of segments in data and cell_imprint are the same
%     if any(strcmp(string(pat_data.segment_id), string(cell_imprint.segment_id)) == 0)
%         [~, id] = ismember(string(cell_imprint.segment_id), string(pat_data.segment_id));
%         % Reorder cell_imprint to match data and metadata
%         cell_imprint = cell_imprint(id,:);
%     end
    
    % Choose atlas
%     for atl = [60, 125]
        % Choose threshold for allowing propagated activity (seconds)
%         for det = [0,1,2]



%% convert channel data to ROI with mapping

[mapping_ch2roi,ch_names,roi_names]=map_to_roi(json_data,atl);


%% detect onset based on channels, then convert to onset ROIs


%[imprint_onsetch,EI_onsetch,PLHG_onsetch]=onset_detection_chbased(cell_imprint,pat_data)
%for sarah


dummy_onset=zeros(75,1);
dummy_onset([1,13:50,55])=1;
imagesc(dummy_onset)

dummy_onset_roi=mapping_ch2roi*dummy_onset
imagesc(dummy_onset_roi)

%% labelled onset convert to ROI
labelled_onset=json_data...
mapping_ch2roi*labelled_onset

%% attach meta data

%%
% 
% %%
%             % Compute onset using imprint, EI, and PLHG
%             onset_output = compute_onset(pat_data, json_data, cell_imprint, sz_count_pat, ...
%                 "atlas",atl, 'det', det, 'reg_thresh', chan_to_roi_thresh);
%         
%             if isempty(onset_output)
%                 fprintf("Patient %s does not meet inclusion criteria \n", patient)
%                 continue
%             elseif isempty(onset_output.Segment_ids{1,1})
%                 continue 
%             end
% 
%             % Compute ROI-based onset for imprint, EI, and PLHG
%             for met = ["imprint", "EI", "PLHG"]
%                 [onset_output] = onset_chan_to_roi(pat_data, json_data, onset_output, ...
%                     'method', met, 'atlas', atl, 'thresh', chan_to_roi_thresh);
%             end
% 
%             onset_output.atl = atl;
%             onset_output.det = sprintf("+%d seconds", det);
%             
%             % Add patient onset_output to final output table
%             if size(onset_output,1) > 0
%                 if exist('final_output', 'var')
%                     final_output = [final_output; onset_output];
%                 else
%                     final_output = onset_output;
%                 end
%             end
% %         end
% 
% 
% %     end
