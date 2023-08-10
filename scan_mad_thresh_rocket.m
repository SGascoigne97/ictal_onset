function scan_mad_thresh_rocket(patient)

data_location = '../Data/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];
% data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
% path_pipeline_exports = [data_location, 'export_ictal'];

% Add paths to all functions required
addpath(genpath('help_functions'))
addpath('lib_biomarkers')
addpath('lib_dataflow')
addpath('sarah_functions')
addpath(genpath('lib'))

% Choose atlas
% Choose threshold for allowing propagated activity (seconds)
min_sz = 1;
min_sz_dur = 7; % minimum seizure duration for inclusion
wind_overlap = 7/8; % Overlap between imprint windows
rec_type = "sec";
rec_thresh = 3;
mad_thresh_scan = 0.2:0.2:20;%20;%0; % Scan from 0 to 20 in steps of 0.2
ons_wind_size = 1;
det = 8;

onset_calc_loc = "imprint_ons"; % Specify folder to store imprint values in
% Need a new folder if using a different subset of the data/different data
% as it will load previous save if folder is not empty
%
% For each patient, compute onset based on imprint, EI, and PLHG
    if exist(sprintf('%s/%s.mat', data_location, patient), 'file')
        load(sprintf('%s/%s.mat', data_location, patient));
        pat_data = data_export;
        clear data_export
        
        % load json_data
        filelist = dir(fullfile(strcat(path_pipeline_exports, "/", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
        filelist = filelist(~[filelist.isdir]);  % remove folders from list
        folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
        load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
         
        % Apply inclusion criteria
        [pat_data, ~, ~] = incl_crit(pat_data, 'min_sz_count', min_sz,...
            'sz_type', ["focal", "sg", "subclin", "N"], 'min_sz_duration', min_sz_dur);
        % Append patient data to data table (if inclusion criteria are met)
        if size(pat_data,1) < min_sz
            fprintf('Patient %s does not meet inclusion criteria \n', patient)
            return
        end
        pat_meta = pat_data(:,1:(end-1));
    
        % Remove seizures from json data that do not meet inclusion criteria
        json_data = json_data(contains(string(extractfield(cat(2,json_data.x_id),...
            'x_oid')), string(pat_data.segment_id)));
        
        % Compute onsets
        mkdir(onset_calc_loc)
    
        % Iterate through MAD thresholds
        pat_tab = table(repmat(patient,size(pat_data, 1),1), string(pat_data.segment_id),...
            'VariableNames', {'patient','segment_id'});
        
    %%
        for mad_thresh = mad_thresh_scan
            % Compute imprint for all seizures
            [pat_data, pat_meta, cell_imprint,  sz_count_pat] = ...
                           calc_imprint_rocket(pat_data, pat_meta, "window_overlap", wind_overlap,...
                                        "folder",onset_calc_loc, "min_sz_count", min_sz,...
                                        "rec_type", rec_type, "rec_thresh",rec_thresh, ...
                                        "mad_thresh", mad_thresh);
            
            %%
            % Check order of segments in data and cell_imprint are the same
            if any(strcmp(string(pat_data.segment_id), string(cell_imprint.segment_id)) == 0)
                [~, id] = ismember(string(cell_imprint.segment_id), string(pat_data.segment_id));
                % Reorder cell_imprint to match data and metadata
                cell_imprint = cell_imprint(id,:);
            end

             %% detect channel-wise onset
            [onset_mat, onset_time] = compute_onset_rocket(pat_data, cell_imprint, 'wdw_sz', ons_wind_size, 'det', det);

            % Create empty array to proportion of windows with activity 
            % detected for each MAD threshold
            pat_ratio = zeros(size(pat_data, 1), 1);
            prop_chan_onset = zeros(size(pat_data, 1), 1);
            
            for sz = 1:size(pat_data, 1)
                % Prop of windows with activity detected
                ratio = sum(sum(cell_imprint{sz,1}{:}))/prod(size(cell_imprint{sz,1}{:}));
                pat_ratio(sz) = ratio;

                % Proportion of channels in onset
                prop_chan_onset(sz) = sum(onset_mat(:,sz))/size(onset_mat,1);

            end
            pat_tab.(sprintf("prop_act_mad_%.1f", mad_thresh)) = pat_ratio;
            pat_tab.(sprintf("prop_chan_ons_mad_%.1f", mad_thresh)) = prop_chan_onset;
            pat_tab.(sprintf("ons_time_mad_%.1f", mad_thresh)) = onset_time';
    
            save(sprintf('output/pat_tab_%s.mat', patient), 'pat_tab')
        
        end
    else
        fprintf('Patient %s does not have saved data \n', patient)
    end

end




% function scan_mad_thresh_rocket(patient)
% 
% data_location = '../Data/'; % replace
% path_pipeline_exports = [data_location, 'export_ictal'];
% % data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
% % path_pipeline_exports = [data_location, 'export_ictal'];
% 
% % Add paths to all functions required
% addpath(genpath('help_functions'))
% addpath('lib_biomarkers')
% addpath('lib_dataflow')
% addpath('sarah_functions')
% addpath(genpath('lib'))
% 
% % Choose atlas
% % Choose threshold for allowing propagated activity (seconds)
% min_sz = 1;
% min_sz_dur = 7; % minimum seizure duration for inclusion
% wind_overlap = 7/8; % Overlap between imprint windows
% rec_type = "sec";
% rec_thresh = 3;
% mad_thresh_scan = 0.2:0.2:20;%0; % Scan from 0 to 20 in steps of 0.5
% 
% onset_calc_loc = "imprint_ons"; % Specify folder to store imprint values in
% % Need a new folder if using a different subset of the data/different data
% % as it will load previous save if folder is not empty
% %
% % For each patient, compute onset based on imprint, EI, and PLHG
%     if exist(sprintf('%s/%s.mat', data_location, patient), 'file')
%         load(sprintf('%s/%s.mat', data_location, patient));
%         pat_data = data_export;
%         clear data_export
%         
%         % load json_data
%         filelist = dir(fullfile(strcat(path_pipeline_exports, "/", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
%         filelist = filelist(~[filelist.isdir]);  % remove folders from list
%         folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
%         load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
%          
% %         % Apply inclusion criteria
% %         [pat_data, ~, ~] = incl_crit(pat_data, 'min_sz_count', min_sz,...
% %             'sz_type', ["focal", "sg", "subclin", "N"], 'min_sz_duration', min_sz_dur);
% %         % Append patient data to data table (if inclusion criteria are met)
% %         if size(pat_data,1) < min_sz
% %             fprintf('Patient %s does not meet inclusion criteria \n', patient)
% %             return
% %         end
%         pat_meta = pat_data(:,1:(end-1));
%     
% %         % Remove seizures from json data that do not meet inclusion criteria
% %         json_data = json_data(contains(string(extractfield(cat(2,json_data.x_id),...
% %             'x_oid')), string(pat_data.segment_id)));
%         
%         % Compute onsets
%         mkdir(onset_calc_loc)
%     
%         % Iterate through MAD thresholds
%         pat_tab = table(repmat(patient,size(pat_data, 1),1), string(pat_data.segment_id),...
%             'VariableNames', {'patient','segment_id'});
%     %%
%         for mad_thresh = mad_thresh_scan
%             % Compute imprint for all seizures
%             [pat_data, pat_meta, cell_imprint,  sz_count_pat] = ...
%                            calc_imprint_rocket(pat_data, pat_meta, "window_overlap", wind_overlap,...
%                                         "folder",onset_calc_loc, "min_sz_count", min_sz,...
%                                         "rec_type", rec_type, "rec_thresh",rec_thresh, "mad_thresh", mad_thresh);
%             %%
%             % Check order of segments in data and cell_imprint are the same
%             if any(strcmp(string(pat_data.segment_id), string(cell_imprint.segment_id)) == 0)
%                 [~, id] = ismember(string(cell_imprint.segment_id), string(pat_data.segment_id));
%                 % Reorder cell_imprint to match data and metadata
%                 cell_imprint = cell_imprint(id,:);
%             end
%             % Create empty array to store activity:no activity ratio for each MAD
%             % threshold
%             pat_ratio = zeros(size(pat_data, 1), 1);
%             
%             for sz = 1:size(pat_data, 1)
%                 % Ratio of ones (activity) against zeros (no activity)
%                 ratio = sum(sum(cell_imprint{sz,1}{:}))/prod(size(cell_imprint{sz,1}{:}));
%                 pat_ratio(sz) = ratio;
%             end
%             pat_tab.(sprintf("mad_%.1f", mad_thresh)) = pat_ratio;
%     
%             save(sprintf('output/pat_tab_%s.mat', patient), 'pat_tab')
%         
%         end
%     else
%         fprintf('Patient %s does not have saved data \n', patient)
%     end
% 
% end
% 
% 
