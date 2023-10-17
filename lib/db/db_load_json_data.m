function json_data = db_load_json_data(file_json,options)
% DB_LOAD_JSON_DATA  Load json export from CNNP iEEG database as a
% structure array in MATLAB. 
%
% In order to form a structure array, the function adds any missing fields
% at the segment level with the value [] (empty array).
%
% By default, channel details of each segments are also converted from a
% structure to a table during the import and all channel details tables are
% standardised to have the same variables.
%
% In order to format segment channel details, channel details MUST be 
% present and contain 2+ channels. Problematic segments will be removed
% from the JSON structure by default.
%
%   data = DB_LOAD_JSON_DATA(file_json) loads a json file, specified by the
%   string or character array file_json, as a structure, json_data.
%
%   data = DB_LOAD_JSON_DATA(...,Name, Value) specifies formatting options
%   for json_data:
%       - 'RemoveNoChannelDetails' (boolean) specifies whether to remove
%       segments that have problematic channel_details (i.e., segments that
%       do not have a channel_details field with 2+ channels). Default:
%       true. 
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
%       - 'Hospital' (cell array) limits the import to segments from
%       hospitals listed in the cell array. Default: {} (no hospital
%       restrictions).
%   
%   For example:
%
%   data = DB_LOAD_JSON_DATA(file_json,'RemoveMissChannelDetails',false,...
%   'ChannelTable',false) will return the JSON structure json_data with all
%   segments and without formatted channel_details fields. 
%
%
%   See also DB_FORMAT_JSON_DATA, DB_RM_SEGM_NO_CHANNEL_DETAILS, 
%   DB_GET_CHANNEL_VARIABLES, DB_IS_FROM_HOSPITAL
%
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022
% January 2023 - added option for limiting segments to certain hospitals
% February 2023 - separated formatting steps into a separate function and
% added option for removing segments whose channel details cannot be
% formatted.
%
% N.B. for future developers: when fields are missing for channels in 
% channel_details, the fuction adds NaNs to numeric arrays (e.g., is_soz) 
% and empty strings to cell arrays. Note this approach could potentially 
% create fields/variables with mixed data types. (e.g., if empty strings 
% are added to ROIids). 
% Any fields/variables with mixed data types must be included in ignore
% list for db_get_channel_variables.

arguments
    file_json char
    options.RemoveNoChannelDetails (1,1) {mustBeNumericOrLogical} = true
    options.ChannelTable (1,1) {mustBeNumericOrLogical} = true
    options.StandardiseChannelTable (1,1) {mustBeNumericOrLogical} = true
    options.Hospital cell = {}
end


% load JSON data
fid = fopen(file_json);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
json_data = jsondecode(str);

% number of segments of data
n_segm = length(json_data);

% check if imported data is a structure ( = all segments have same fields)
if isstruct(json_data)
    disp('JSON data imported as structure.')
elseif iscell(json_data)
    disp('JSON data imported as cell array of structures - adding missing fields and converting to structure.')
    
    % all structures must have same fields in the same order in order to use
    % cell2mat() to combine as array
    
    % get complete list of fields of each stuct
    json_fields = {};
    for i=1:n_segm
        json_fields = unique([json_fields; fieldnames(json_data{i})],'stable');
    end
    
    % add any missing fields and ensure same field order
    for i=1:n_segm
        % find missing fields in each segment
        fields_i = fieldnames(json_data{i});            % segment i's fields
        miss_fields = setdiff(json_fields,fields_i);    % missing fields in segment i
        n_miss = length(miss_fields);                   % n missing
        
        % add any missing fields (as empty array)
        if n_miss > 0
            for j=1:n_miss
                json_data{i}.(miss_fields{j}) = [];
            end
        end
        
        % ensure field order is the same as in json_fields
        json_data{i} = orderfields(json_data{i},json_fields);
    end
    
    % now ready to turn into a structure array
    json_data = cell2mat(json_data);
else
    error('JSON data imported in unexpected format - import code will need to be modified to handle new format.')
end

% restrict import to segments from specific hospitals if requested
if ~isempty(options.Hospital)
    disp('RESTRICTING TO SEGMENTS FROM SPECIFIED HOSPITALS')
    keep_bool = db_is_from_hospital(json_data,options.Hospital);
    disp(['Removing ' num2str(sum(keep_bool==0)) ' segments from unrequested hospitals']) % number of segments removed
    json_data = json_data(keep_bool);                       % subset json_data
    n_segm = length(json_data);                             % new number of segments
end

%%%%% REMOVING SEGMENTS WITH INSUFFICIENT CHANNEL DETAILS %%%%%
% = segments that cannot be formatted due to 0 or 1 channels in channel_details
if options.RemoveNoChannelDetails
disp('REMOVING SEGMENTS WITH 0 OR 1 CHANNELS IN CHANNEL DETAILS')

    % remove
    [json_data, rm_bool] = db_rm_segm_no_channel_details(json_data);
    
    % number of segments removed
    disp(['Removed ' num2str(sum(rm_bool)) ' segments with insufficient channel details'])
end

%%%%% FORMATTING CHANNEL DETAILS USING DB_FORMAT_JSON_DATA %%%%%
json_data = db_format_json_data(json_data,...
    'ChannelTable',options.ChannelTable,...
    'StandardiseChannelTable',options.StandardiseChannelTable);

