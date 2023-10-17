function json_data = pre_filt_chan_by_feature(json_data,keep_values)
% PRE_FILT_CHAN_BY_FEATURE Add variable to JSON structure indicating which
% channels meet inclusion criteria. Can filter channels based on variables
% in channel details tables.
%
%   json_data = PRE_FILT_CHAN_BY_FEATURE(json_data, keep_values) modifies
%   json_data by 
%       1) adding an additional logical variable to each channel 
%       details table, keep_chan, indicating which channels match the inclusion
%       criteria specified by keep_values. keep_values is a structure with
%       fields matching the variables in the channel details tables, and each
%       field contains the values permitted for that variable.
%       and
%       2) listing names of channels to remove in an additional field,
%       bad_chan_filtered, in each segment (json_data.bad_chan_filtered)
%
%   EXAMPLES:
%
%   EXAMPLE 1: Keep spared depth electrodes; must also be outside of the
%   seizure onset zone if SOZ info known
%       keep_values = struct();
%       keep_values.contact_type = 'depth'; % keep depth electrodes
%       keep_values.is_resected5 = 0; % keep spared channels (based on 5mm distance threshold to mask)
%       keep_values.is_soz = [0,NaN]; % keep channels if confirmed not SOZ (0) or unknown if SOZ (NaN)
%
%   EXAMPLE 2: Keep grid and strip electrodes; do not filter based on any
%   other criteria
%       keep_values = struct();
%       keep_values.contact_type = {'grid','strip'}; % keep grid and strip electrodes
%
% Channels must meet ALL inclusion criteria in keep_values. 
%
% A warning will be thrown if keep_values contains any fields that are not
% variables in the channel details tables, and those variables will not be
% used for restricting channels. User will be asked to input y to continue
% despite missing variables.
%
% Recommend using DB_GET_CHANNEL_DETAILS_VALUES to determine which values
% are present in channel details in your data. 
%
% Note that this function is responsible for marking any scalp EEG
% electrodes for removing (using contact_type variable). 
%
% See also DB_GET_CHANNEL_VARIABLES, DB_GET_CHANNEL_DETAILS
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
    keep_values struct
end

disp('MARKING CHANNELS TO REMOVE BASED ON CHANNEL FEATURES')

% variables to use for channel inclusion criteria
crit_var = fieldnames(keep_values); % criteria variables
n_crit_var = length(crit_var);      % number of variables

% check that all keep_values fields are also variables in channel details tables
% throws warning if mismatch (protects against mistakes due to, e.g., typos)
all_var = db_get_channel_variables(json_data);
is_match = ismember(crit_var, all_var); 
if sum(is_match) < n_crit_var
    miss_var = cellfun(@string,crit_var(is_match == 0)) + ' ';
    warning(['Variables not in channel details ( ' miss_var{:} ') will not be used for channel inclusion criteria'])
    txt = input('Do you wish to continue? Type y or n and press Enter',"s");
    switch txt
        case 'n'
            error('Function stopped due to missing variable at request of user')
        case 'y'
            crit_var = crit_var(is_match);
            n_crit_var = length(crit_var);
            disp('Missing variables will not be used.')
        otherwise
            error('input must be y or n')
    end
end

% mark which channels match criteria
n_segm = length(json_data);
for i=1:n_segm
    chan_tbl = json_data(i).channel_details;    % segment's channel details table
    n_chan = height(chan_tbl);                  %  number of channels
    
    % initialise array for storing whether channel is a match for each criteria
    match_crit = zeros(n_chan,n_crit_var); 
    
    % check each criteria
    for j=1:n_crit_var
        chan_crit = keep_values.(crit_var{j});  % allowed values for variable
        chan_data = chan_tbl.(crit_var{j});     % channel data for variable 
        match_crit(:,j) = ismember(chan_data,chan_crit); 
        
        % ismember treats NaN values as distinct - need to separately check
        % for NaN values if NaNs permitted
        if isnumeric(chan_crit)             % if numeric, check for nans
            if sum(isnan(chan_crit))>0      % check if nans allowed
                is_nan = isnan(chan_data);  % nans in channel data
                match_crit(:,j) = match_crit(:,j) | is_nan; % add nans to matching vector
            end
        end
        
    end
    
    % add variable for whether each channel meets inclusion criteria
    json_data(i).channel_details.keep_chan = sum(match_crit,2) == n_crit_var;
    
    % display number of channels kept in analysis
    n_kept = sum(json_data(i).channel_details.keep_chan);
    disp(['segment ' num2str(i) ': ' num2str(n_kept) '/' num2str(n_chan) ' channels in channel details have requested features'])
    
    % list names of channels to remove and save in json structure
    rm_chan = json_data(i).channel_details.chan_name(json_data(i).channel_details.keep_chan == 0);
    json_data(i).bad_chan_filtered = rm_chan;
    
end