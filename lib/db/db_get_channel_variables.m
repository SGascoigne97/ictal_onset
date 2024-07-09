function channel_variables = db_get_channel_variables(json_data)
% DB_GET_CHANNEL_VARIABLES Get list of all variables included in the
% channels details tables in the JSON structure.
%
%   channel_variables = DB_GET_CHANNEL_VARIABLES(json_data) returns a cell
%   array, channel_variables, of all variables found in the channel details
%   tables in the JSON structure json_data. Variables are included even if
%   they are not found in all channel details tables.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
end

% number of iEEG segments
n_segm = length(json_data);

% initialise array
channel_variables = {};

for i=1:n_segm
    
    % variables in channels details table
    chan_var = json_data(i).channel_details.Properties.VariableNames;
    
    % add to list of all variables; only keep unique
    channel_variables = unique([channel_variables chan_var],'stable');
end