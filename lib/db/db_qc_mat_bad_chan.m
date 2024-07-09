function json_data = db_qc_mat_bad_chan(json_data,mat_path)
% DB_QC_MAT_BAD_CHAN Make sure that eeg_bad_chan variables with no entries
% are empty cell arrays and that all eeg_bad_chan and eeg_channel variables
% are column vectors.
%
%   DB_QC_MAT_BAD_CHAN(json_data,mat_path) formats all eeg_bad_chan and
%   eeg_channels variables of the iEEG segments in the JSON structure 
%   json_data. The iEEG workspaces are stored in the folder mat_path. 
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022
% modified February 2023 to also format eeg_channels - may want to revise
% name in future to indicate more general usage.

arguments
    json_data struct
    mat_path char
end

% number of iEEG segments
n_segm = length(json_data);

% get file names of all iEEG segments
json_eeg_fn = db_all_json_eeg_fn(json_data);

disp('Checking and converting format of eeg_bad_chan...')

% Quality control of formatting of segments with no bad channels
for i=1:n_segm
    clearvars eeg_bad_chan eeg_channels
    % load
    load([mat_path '/' json_eeg_fn{i}],'eeg_bad_chan','eeg_channels','eeg_fs');

    if ~exist('eeg_bad_chan', 'var')
        eeg_bad_chan = {};
    end
   % make empty cell array if 1x1 cell array of 0x0 char
    if length(eeg_bad_chan)==1
        if isempty(eeg_bad_chan{1})
            eeg_bad_chan = {};
            disp(json_eeg_fn{i})
        end
    end
    
    % also ensure output is column vector
    if size(eeg_bad_chan,2) > 1
        eeg_bad_chan = eeg_bad_chan';
    end

    if ~exist('eeg_channels', 'var')
        eeg_channels = {};
    end
    
    % make sure eeg_channels is a column vector
    if size(eeg_channels,2)>1
        eeg_channels = eeg_channels';
    end
    
    % add to json_data
    json_data(i).eeg_bad_chan = eeg_bad_chan;
    json_data(i).eeg_channels = eeg_channels;
    json_data(i).eeg_fs = eeg_fs;
end