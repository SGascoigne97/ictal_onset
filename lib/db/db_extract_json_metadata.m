function json_metadata = db_extract_json_metadata(json_data,metadata_type)
% DB_EXTRACT_JSON_METADATA Extract one type of metadata for all iEEG
% segments in the JSON structure.
%
%   json_metdata = DB_EXTRACT_JSON_METADATA(json_data,metadata_type)
%   extracts a list of metadata, json_metdata, for all segments in the JSON
%   structure json_data for the specified metadata type, metadata_type. 
%   json_metdata is a cell array (length equal to the number of segments)
%   unless otherwise specified below.
%
%   Options for metadata type:
%       'x_oid': segment JSON ID
%       'segm_num': segment number (see db_all_segm_num) (numeric)
%       'eeg_duration': length of segment in seconds (numeric)
%       'Hospital': hospital where data was collected
%       'patID': patient ID
%       'patientDOB': patient date of birth (year) (numeric)
%       'patientSex': patient sex (M/F)
%       'patientAge_onset': patient age of onset of epilepsy (numeric)
%       'patientAge_exam_est': estimated patient age at exam (numeric)
%           (derived from exam_date and patientDOB)
%       'exam_id': database ID of exam that segment belongs to
%       'exam_date': date of iEEG recording
%       'outcome_ILAE1': first ILAE outcome listed for patient (numeric)
%       'outcome_ILAE': all ILAE outcomes listed for patient
%       'outcome_Engel1': first Engel outcome listed for patient (numeric)
%       'outcome_Engel': all Engel outcomes listed for patient
%       'outcome_year': corresponding years of ILAE or Engel outcomes
%       'outcome_date': corresponding dates of ILAE or Engel outcomes
%       'outcome1_minus_exam_year': first outcome year minus exam year
%           (derived from outcome_year and exam_date)
%       'op_type': type of operation performed
%       'op_lobe': lobe(s) operaton was performed on (derived from op_type)
%       'op_side': hemisphere of operation
%       'op_pathology': pathology label from operation
%
%   See also DB_EXTRACT_JSON_METADATA_TABLE, DB_EXTRACT_AND_SAVE_JSON_METADATA,
%   DB_ALL_JSON_ID, DB_ALL_JSON_HOSPITAL, DB_ALL_SEGM_NUM
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023
%
% Development note: if add additional metadata types, please also add to
% default variables for db_extract_json_metadata_table

arguments
    json_data struct
    metadata_type char {mustBeMember(metadata_type,{...
        'x_oid','segm_num',...
        'eeg_duration','Hospital','patID','patientDOB',...
        'patientSex','patientAge_onset','patientAge_exam_est',...
        'exam_id','exam_date','outcome_ILAE1','outcome_ILAE',...
        'outcome_Engel1','outcome_Engel','outcome_year','outcome_date',...
        'outcome1_minus_exam_year',...
        'op_type','op_lobe','op_side','op_pathology'})}
end
    
% number of iEEG segments
n_segm = length(json_data);

% initialise cell array
json_metadata = cell(n_segm,1);

% get requested metadata
% I've noted which types are assumed to be present for all segments and
% which types had pre-existing functions to extract that data.
switch metadata_type
    case 'x_oid'        % assumed present
        json_metadata = db_all_json_id(json_data);          % pre-existing function
        
    case 'segm_num'     % assumed present
        json_metadata = db_all_segm_num(json_data);  % pre-existing function
        
    case 'eeg_duration' % assumed present
        for i=1:n_segm
            json_metadata{i} = json_data(i).eeg_duration;
        end
        
    case 'Hospital'     % assumed present
        json_metadata = db_all_json_hospital(json_data) ;   % pre-existing function
        
    case 'patID'        % assumed present
        for i=1:n_segm
            json_metadata{i} = json_data(i).patient_details.patID;
        end
        
    case 'patientDOB'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).patient_details.DOB;
            catch
                json_metadata{i} = NaN; % NaN if missing since converted to numeric
            end
        end
        
    case 'patientSex'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).patient_details.Sex;
            catch
                json_metadata{i} = '';
            end
        end
        

    case 'patientAge_onset'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).patient_details.Age_onset;
            catch
                json_metadata{i} = NaN; % NaN if missing since converted to numeric
            end
        end
        
    case 'patientAge_exam_est' % derived from other metadata
        dob_year = db_extract_json_metadata(json_data,'patientDOB');                % all DOB year
        exam_date = db_extract_json_metadata(json_data,'exam_date');                % all exam date
        exam_year = cellfun(@(x) x(1:4),exam_date,'UniformOutput',false);           % exam year
        exam_year = cell2mat(cellfun(@str2num,exam_year,'UniformOutput',false));    % as numeric
        
        json_metadata = exam_year - dob_year; % difference between exam year and DOB year
        
    case 'exam_id' % assumed present
        for i=1:n_segm 
            json_metadata{i} = json_data(i).exam_id.x_oid;
        end
        
    case 'exam_date'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).exam_details.exam_date.x_date;
            catch
                json_metadata{i} = '';
            end
        end
        
    case 'outcome_ILAE1'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.outcome_ILAE(1);
            catch
                json_metadata{i} = NaN; % NaN if missing since converted to numeric
            end
        end
        
    case 'outcome_ILAE'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.outcome_ILAE;
            catch
                json_metadata{i} = []; % empty array since data is an array
            end
        end
        
     case 'outcome_Engel1'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.outcome_Engel(1);
            catch
                json_metadata{i} = NaN; % NaN if missing since converted to numeric
            end
        end   
               
    case 'outcome_Engel'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.outcome_Engel;
            catch
                json_metadata{i} = []; % empty array since data is an array
            end
        end   
        
    case 'outcome_year'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.outcome_year;
            catch
                json_metadata{i} = []; % empty array since data is an array
            end
        end
        
    case 'outcome_date'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.outcome_date;
            catch
                json_metadata{i} = []; % empty array since data is an array
            end
        end
        
    case 'outcome1_minus_exam_year' % derived from other metadata
        % exam year
        exam_date = db_extract_json_metadata(json_data,'exam_date');                % all exam date
        exam_year = cellfun(@(x) x(1:4),exam_date,'UniformOutput',false);           % exam year
        exam_year = cell2mat(cellfun(@str2num,exam_year,'UniformOutput',false));    % as numeric
        
        % outcome year, first only 
        outcome_year = zeros(n_segm,1);
        for i=1:n_segm
            try
                outcome_year(i) = json_data(i).treatment_details.outcome_year(1);
            catch
                outcome_year(i) = NaN;
            end
        end            
        json_metadata = outcome_year - exam_year; % difference between first outcome year and exam year  
        
    case 'op_type'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.op_type;
            catch
                json_metadata{i} = '';
            end
        end
        
    case 'op_lobe' % derived from other metadata
        op_type = db_extract_json_metadata(json_data,'op_type');
        % remove non-lobe abbreviations and extra spaces
        op_lobe = strrep(op_type,'Lx','');
        op_lobe = strrep(op_lobe,'Lesx','');
        op_lobe = strip(strrep(op_lobe,' ',''));
        
        json_metadata = op_lobe;
        
    case 'op_side'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.op_side;
            catch
                json_metadata{i} = '';
            end
        end
        
    case 'op_pathology'
        for i=1:n_segm
            try
                json_metadata{i} = json_data(i).treatment_details.op_pathology;
            catch
                json_metadata{i} = '';
            end
        end 
        
end

% convert certain metadata to numeric array
% segm_num, patientAge_exam_est, and outcome1_minus_exam_year already numeric - no need to convert
switch metadata_type
    case {'eeg_duration','patientAge_onset',...
            'patientDOB','outcome_ILAE1',...
            'outcome_Engel1'}
        json_metadata = cell2mat(json_metadata);
end

