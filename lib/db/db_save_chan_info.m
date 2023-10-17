%% Helper function for extracting channel information for each patient
function [pat_chan_details] = db_save_chan_info(patient, json_data, save_location, opts)
% Given patient ID and location of data, this function will extract and
% store channel information for this patient
% 
% input:
%   - patient: patient ID (without UCLH)
%   - data_location: filepath to where patient json data is stored

% output
%   - pat_chan_details: table of channel details only for channels
%   remaining in the data after preprocessing

    arguments
        patient (1,1) string % Patient ID
        json_data % Structure containing patient data
      %  data_location (1,1) string % Path to where patient data is stored
        save_location (1,1) string % Path to where you want channels tables to be saved
        opts.save_tab (1,1) double {mustBeNumericOrLogical} = 1 % If 1, save table
    end
    
    % Add in optional arguments
    save_tab = opts.save_tab;

%     % Load json data
%     filelist = dir(fullfile(strcat(data_location, "/UCLH", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
%     filelist = filelist(~[filelist.isdir]);  % remove folders from list
%     folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
%     load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
    
    % Pull out channel information only for channels included after
    % preprocessing

    if isfield(json_data(1), 'pre_eeg_channels') 
        pat_chan_details = json_data(1).channel_details(find(...
            ismember(strrep(json_data(1).channel_details.chan_name, ' ',''),...
            strrep(json_data(1).pre_eeg_channels, ' ',''))),:);
        if save_tab == 1
            mkdir(sprintf('%schannels', save_location))
            save(sprintf("%schannels/channels_%s.mat",save_location, patient), 'pat_chan_details', '-v7.3');
        end
    else
        fprintf("Patient %s does not have pre_eeg_channels column \n", patient)
    end
 
end
