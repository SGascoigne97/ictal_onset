function [segm_num_new,change_bool] = db_renumber_segm_num(json_data,segm_num)
% DB_RENUMBER_SEGM_NUM Ensures segments of each exam are numbered from 1 to
% n (n = number of segments); i.e., no numbers are skipped.
%
% This function is helpful if you want to ensure that you are analysing the
% largest possible number of exams when subsetting by segm_num == 1 (and
% so on); e.g., you will not miss an exam that only has one segment
% labelled with segment number 2. 
%
%   [segm_num_new,change_bool] = DB_RENUMBER_SEGM_NUM(json_data,segm_num)
%   takes a vector of the current segment numbers, segm_num, along with the
%   associated JSON structure json_data, and returns a new vector of
%   segment numbers, segm_num_new, with each exam's segments numbered from
%   1 to n segments. It also returns change_bool, a boolean of whether a
%   segment's number changed.
%
%   Segments are renumbered using their database ID alphanumeric order.
%
%   See also DB_ALL_SEGM_NUM
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    json_data struct
    segm_num (:,1) {mustBeNumeric}
end

% check that size of json_data matches size of segm_num
if length(segm_num) ~= length(json_data)
    error('Sizes of segm_num and json_data do not match')
end

% get array of correspoding exams that segments belong to
json_exam_id = db_extract_json_metadata(json_data,'exam_id');

% also get segment identifiers - will use to renumber segments if necessary
json_segm_id = db_extract_json_metadata(json_data,'x_oid');

% list of unique exam ids 
unique_exam = unique(json_exam_id);
n_exam = length(unique_exam);

% initialise segm_num_new with current numbers
segm_num_new = segm_num;

% loop through exams and check that segment numbers are 1:n
for i=1:n_exam
    
    % exam id
    eid = unique_exam{i};
    
    % exam segments and the corresponding segment num
    exam_segm_idx = find(ismember(json_exam_id,eid));
    exam_segm_num = segm_num(exam_segm_idx);
    n_exam_segm = length(exam_segm_idx);
    
    % check that numbers are 1:n
    if ~isequal(sort(exam_segm_num)',1:n_exam_segm)
        disp(['renumbering segments of exam ' eid])
        
        % get segment ids
        sid = json_segm_id(exam_segm_idx);
        
        % get ranking of segment ids
        [~,~,s_rank] = unique(sid);
        
        % use ranks as segment numbers
        segm_num_new(exam_segm_idx) = s_rank;
        
    end
    
end

% which segment numbers were changed
change_bool = segm_num_new ~= segm_num;
disp(['Changed numbers of ' num2str(sum(change_bool)) ' segments'])

end

        
    
