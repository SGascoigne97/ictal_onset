function segm_num = db_all_segm_num(json_data)
% DB_ALL_SEGM_NUM Extract/assigns segment numbers from the original
% iEEG workspace names.
%
%   segm_num = DB_ALL_SEGM_NUM(json_data) extracts or assigns segment
%   numbers (segm_num) for all iEEG workspace names in the JSON structure 
%   json_data. Segment numbers are preferentially pulled from the text 
%   between "segm" and ".mat" in each workspace name and are returned as a 
%   numeric. If the filename does not include segm#.mat, then segment
%   numbers are instead assigned by determining the number of segments in
%   each exam and assigning numbers in the alphanumeric order of the
%   segment database IDs. 
%
%   See also DB_EXTRACT_JSON_METADATA, DB_RENUMBER_SEGM_NUM
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023 
% modified February 2023

arguments
    json_data struct
end

% get segment filenames
json_fn = db_all_json_eeg_fn(json_data);

% number of iEEG segments
n_segm = length(json_fn);

% initialise arrays
segm_num = nan(n_segm,1);

% locations of eeg/, segm and .mat in each filename
segm_loc = strfind(json_fn,'segm');
mat_loc = strfind(json_fn,'.mat');


% first pass: extract text between "segm" and ".mat" in each filename
% if no "segm" or segm not followed by number, then will leave NaN for now

for i=1:n_segm
    fn = json_fn{i};
    if ~isempty(segm_loc{i})
        sloc_i = segm_loc{i}(end); % end ensures last occurrence if multiple
        mloc_i = mat_loc{i}(end);
        segm_num(i) = str2double(fn((sloc_i+4):(mloc_i-1))); % expected to be numeric - will return NaN if letter
    end
end


% second pass: assign segment numbers to segments without segment labels
% will base number on segment database id

% get array of correspoding exams that segments belong to
json_exam_id = db_extract_json_metadata(json_data,'exam_id');

% also get segment identifiers - will use to number segments if not already
% numbered
json_segm_id = db_extract_json_metadata(json_data,'x_oid');

% list of unique exam ids that have nan segment numbers
unique_exam = unique(json_exam_id(isnan(segm_num)));
n_exam = length(unique_exam);

% loop through exams 
for i=1:n_exam
    
    % exam id
    eid = unique_exam{i};
    
    % exam segments and the corresponding segment ids
    exam_segm_idx = find(ismember(json_exam_id,eid));
    sid = json_segm_id(exam_segm_idx);
    
    % get ranking of segment ids
    [~,~,s_rank] = unique(sid);
    
    % use ranks as segment numbers 
    segm_num(exam_segm_idx) = s_rank;
    
end


% check that no overlap between segment numbers in same exam
% get array of correspoding exams that segments belong to

% list of all unique exam ids 
unique_exam = unique(json_exam_id);
n_exam = length(unique_exam);

% loop through exams and check segment numbers
for i=1:n_exam
    
    % exam id
    eid = unique_exam{i};
    
    % exam segments and the corresponding segment ids
    exam_segm_idx = find(ismember(json_exam_id,eid));
    sid = json_segm_id(exam_segm_idx);
    
    % unique segment numbers 
    snum = unique(segm_num(exam_segm_idx));
    
    % check that length matches
    if length(sid) ~= length(snum)
        error(['Segments in exam ' eid ' does not have unique segment numbers'])
    end    
end

    


