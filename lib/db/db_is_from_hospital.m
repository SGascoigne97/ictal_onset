function is_from_hospital = db_is_from_hospital(json_data,hospital)
% DB_IS_FROM_HOSPITAL Returns a logical vector of which segment in the JSON
% structure are from the requested hospital(s).
%
% hospital = 'RAM' is a special case that returns all RAM sites by
% identifying segments whose hospitals are RAM + another letter. 
%
%   is_from_hospital = db_is_from_hospital(json_data,hospital) returns a
%   logical column vector, is_from_hospital, of the segments in the JSON
%   structure json_data that are from the hospital(s) specified by
%   hospital. The variable hospital can be a character array of a single
%   hospital or a cell array of character arrays of multiple hospitals. 
%
%   If hospital contains 'RAM', segments from all RAM hospital sites
%   (which are all RAM plus another letter) are identified.
%
%   See also DB_ALL_JSON_HOSPITAL
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023 

arguments
    json_data struct
    hospital 
end

% check hospital is a cell or character array
if ~iscell(hospital) && ~ischar(hospital)
    error('hospital must be a cell array of multiple hospitals or character array of one hospital')
    
end

% all hospitals in json structure
json_hospital = db_all_json_hospital(json_data);  

% exact matches to hospital(s)
is_from_hospital = ismember(json_hospital,hospital);   

% if RAM resquested, also check for any RAM sites
if ismember({'RAM'},hospital)
    
    % match
    is_ram = regexp(json_hospital,'RAM[A-Z]');
    is_ram = cellfun(@(x) isequal(x,1),is_ram);
    
    % add to is_from_hospital
    is_from_hospital = (is_from_hospital + is_ram)>=1;
    
end