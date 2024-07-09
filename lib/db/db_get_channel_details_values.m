function [val_cell,val_num,val_struct] = db_get_channel_details_values(json_data,options)
% DB_GET_CHANNEL_DETAILS_VALUES Get all values for channel details varables
% that are in a json data structure.
%
% Values are stored in three variables based on the variable data type.
%
%   [val_cell,val_num,val_struct] = 
%   DB_GET_CHANNEL_DETAILS_VALUES(json_data) returns three variables 
%   containing variable values in json_data - one for variables containing
%   cell arrays ("val_cell"), one for variables containing numeric/logical
%   arrays ("val_num"), and one for variables containing structures
%   ("val_struct", returns structure field names only). Each variable is
%   formatted as a structure, with the field names equal to the variable
%   names in the channel details tables. Fields are empty if that variable
%   does not contain the corresponding data type; e.g., numeric array
%   values will only appear in val_num.
%
%   [val_cell,val_num,val_struct] = 
%   DB_GET_CHANNEL_DETAILS_VALUES(json_data,'IgnoreVar',IgnoreVar) excludes
%   the variable names in the cell array IgnoreVar from the exported value
%   variables. By default, IgnoreVar = {'x_id','patient_id','exam_id',
%   'chan_name','location_orig','ROIids','ROIname'}. Note that this 
%   function does not work on variables that are cell arrays of arrays;
%   thus, Ignore Var should include {'elecs_pial','location_orig','ROIids',
%   'ROIname'} at a minimum.
%
% See also DB_GET_CHANNEL_VARIABLES
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
    options.IgnoreVar = {'x_id','patient_id','exam_id','chan_name',...
        'location_orig','location_tal','location_mni','elecs_pial','ROIids','ROIname'};
end

% number of segments
n_segm = size(json_data,1);

% list of all variables in channel details
all_var = db_get_channel_variables(json_data);

% remove variables in ignore list
all_var = setdiff(all_var,options.IgnoreVar,'stable');

% get all values for each variables
n_var = length(all_var);    % number of variables
val_cell = struct();        % for storing values from cell arrays
val_num = struct();         % for storing values in numeric/logical arrays
val_struct = struct();      % for storing field names of structures

% initialise empty array for each variable
for j=1:n_var
    val_cell.(all_var{j}) = {};
    val_num.(all_var{j}) = [];
    val_struct.(all_var{j}) = {};
end

% get values from each segment and keep unique values
for i=1:n_segm
    % segment channel details
    chan_tbl = json_data(i).channel_details;
    
    % values for each variable in that segment 
    for j=1:n_var
        val_ij = chan_tbl.(all_var{j});
        
        try
            % structure fields
            if isstruct(val_ij)
                fields_ij = fieldnames(val_ij);
                val_struct.(all_var{j}) = ...
                    unique([val_struct.(all_var{j}); fields_ij]);
                
                % cell unique values
            elseif iscell(val_ij)
                
                val_cell.(all_var{j}) = ...
                    unique([val_cell.(all_var{j}); val_ij]);
                
                % numeric/logical array unique values
            else
                val_num.(all_var{j}) = ...
                    unique([val_num.(all_var{j}); val_ij]);
            end
        catch
            error(['Cannot get values for variables for ' all_var{j} ' - may need to add variable name to IgnoreVar array when running function'])
        end
    end
end

% remove any duplicate NaNs from numeric arrays 
for j=1:n_var
    val_num.(all_var{j}) = rm_duplicate_nan(val_num.(all_var{j}));
end


function x = rm_duplicate_nan(x)
% remove duplicate nans after calling function "unique" (treats nans as
% distinct)

% check for nans
if sum(isnan(x)) > 0
    % remove all
    x(isnan(x)) = [];
    % add nan to end
    x(end+1) = NaN;
end