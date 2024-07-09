function json_data = db_format_json_data(json_data,options)
% DB_FORMAT_JSON_DATA Formats the channel details of imported JSON data.
%
% This function is called in db_load_json_data by default, but it can also
% be applied to JSON structures that were imported without formatting
% channel details. 
%
% See example code below.
%
%   json_data = DB_FORMAT_JSON_DATA(json_data) formats the JSON structure
%   json_data so that 1) each segment's channel_details is a table, and 2)
%   each channel_details table contains the same variables.
% 
%   The formatting settings can be modified using the following name/value
%   pair arguments:
%       - 'ChannelTable' (boolean) specifies whether to convert each
%       channel_details entry to a table (default: true). All channels are
%       forced to have the same variables in this conversion; any missing
%       variables for a channel are given NaN values.
%       - 'StandardiseChannelTable' (boolean) specifies whether to
%       standardise the channel_details channel tables such that each
%       table has the same variables across all segments in json_data. Any
%       missing variables are given empty strings if a cell array and NaN
%       values otherwise. 'ChannelTable' must also be set to true to 
%       perform this step. Default: true 
%  
%
%   See also DB_LOAD_JSON_DATA, DB_GET_CHANNEL_VARIABLES, 
%   DB_RM_SEGM_NO_CHANNEL_DETAILS
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    json_data struct
    options.ChannelTable (1,1) {mustBeNumericOrLogical} = true
    options.StandardiseChannelTable (1,1) {mustBeNumericOrLogical} = true
end

% do not standardise if ChannelTable = false
if options.StandardiseChannelTable
    if ~options.ChannelTable
        warning('Will not standardise channel details; must first be converted to tables (ChannelTable = true)')
        options.StandardiseChannelTable = false;
    end
end

% number of segments of data
n_segm = length(json_data);

% convert channel details to tables if requested
if options.ChannelTable
    disp('CONVERTING CHANNEL DETAILS TO TABLES')
    for i=1:n_segm
        
        % segment channel details
        chan_det = json_data(i).channel_details;
        
        % convert structure to table
        if isstruct(chan_det)
            
            json_data(i).channel_details = struct2table(chan_det);
            
        else  % stored as cell = not all channels have same fields
            
            % convert each channel's structure to table
            chan_det = cellfun(@db_chan_struct2tbl,chan_det,'UniformOutput',false);
            
            % number of channels
            n_chan = length(chan_det);
            
            % start table with first channel's table
            tbl_all = chan_det{1};
            
            % keep track of missing variables to report
            miss_var = {};
            
            % add each additional channel to table
            for j=2:n_chan
                
                % channel i's table
                tbl_i = chan_det{j};
                
                % find variables present in one table but not the other
                all_miss = setdiff(tbl_i.Properties.VariableNames,tbl_all.Properties.VariableNames);
                i_miss = setdiff(tbl_all.Properties.VariableNames,tbl_i.Properties.VariableNames);
                
                % add nans for any missing variables
                if ~isempty(all_miss)
                    tbl_all = [tbl_all array2table(nan(height(tbl_all),numel(all_miss)),'VariableNames',all_miss)];
                    miss_var = unique([miss_var all_miss]);
                end
                if ~isempty(i_miss)
                    tbl_i = [tbl_i array2table(nan(height(tbl_i),numel(i_miss)),'VariableNames',i_miss)];
                    miss_var = unique([miss_var i_miss]);
                end
                
                % concatenate
                tbl_all = [tbl_all; tbl_i];
                
            end
            
            % convert any NaN entries in cell arrays to empty strings 
            % (otherwise will cause type errors later on)
            if ~isempty(miss_var) 
                for j=1:length(miss_var)
                    if iscell(tbl_all.(miss_var{j}))
                        nan_bool = cellfun(@any_nan,tbl_all.(miss_var{j}));
                        tbl_all.(miss_var{j})(nan_bool) = repmat({''},sum(nan_bool),1);
                    end
                end
            end
            
            % print added variables
            miss_var = cellfun(@string,miss_var) + ' ';
            disp([num2str(i) ': variables added to some channels: ' miss_var{:}])
            
            % replace in json_data
            json_data(i).channel_details = tbl_all;
        end
        
    end
end


if options.StandardiseChannelTable
    
    disp('STANDARDISING CHANNEL TABLE VARIABLES ACROSS ALL SEGMENTS')
    
    % first get complete list of possible channel variables
    channel_variables = db_get_channel_variables(json_data);
    
    % determine whether each variable is stored as a cell array (in any
    % table)
    n_var = length(channel_variables);
    is_cell = zeros(1,n_var);
    
    for i=1:n_segm
        chan_tbl = json_data(i).channel_details;
        for j=1:n_var
            if ismember(channel_variables{j},chan_tbl.Properties.VariableNames)
                if iscell(chan_tbl.(channel_variables{j}))
                    is_cell(j) = 1;
                end
            end
        end
    end  
        
    % standardise
    for i=1:n_segm
        
        % variables in channels details table
        chan_tbl = json_data(i).channel_details;
        chan_var = chan_tbl.Properties.VariableNames;
        
        % variables in complete list that are missing
        [miss_var,idx] = setdiff(channel_variables,chan_var,'stable');
        miss_is_cell = is_cell(idx);
        
        % add nans or empty strings for any missing variables
        if ~isempty(miss_var)
            n_miss = length(miss_var);
            for j=1:n_miss
                if miss_is_cell(j) == 1
                    chan_tbl = [chan_tbl array2table(repmat({''},height(chan_tbl),1),...
                        'VariableNames',miss_var(j))];
                else
                    chan_tbl = [chan_tbl array2table(nan(height(chan_tbl),1),...
                        'VariableNames',miss_var(j))];
                end
            end
            
            % put new table in json_data
            json_data(i).channel_details = chan_tbl;
            
            % display
            miss_var = cellfun(@string,miss_var) + ' ';
            disp([num2str(i) ': added NaNs or empty strings for ' miss_var{:}])
            
        else
            disp([num2str(i) ': all variables present'])
        end
        
    end
end



% convert each channel's structure to a table
function chan_tbl = db_chan_struct2tbl(chan_struct)

chan_tbl = struct2table(chan_struct,'AsArray',true);


function has_nan = any_nan(x)

has_nan = sum(isnan(x)) > 0;