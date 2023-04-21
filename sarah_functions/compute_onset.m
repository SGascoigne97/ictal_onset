function [onset_output] = compute_onset(data_tbl, json_data, cell_imprint, sz_count_pat)
% Relabel patients to match Seizure Severity paper (Gascoigne et al., 2023)

% At the moment, the code only works for one or two seizure type(s) at a 
% time, this can be adjusted to inmclude multiple seizure types 

% input:
%   - data_tbl: full data table
%   - cell_imprint: table of imprints for all seizures included in data table

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl
        json_data
        cell_imprint 
        sz_count_pat
    end
    
    % We could include an optional argument to determine which onset detection
    % methods we will use
    % Additional arguents to change parameters in detection methods (e.g.,
    % lambda and v in EI)
    
    % Create table to store output
    colNames = {'Patient_id', 'Segment_ids', 'channel_names', 'ROI_ids',...
        'labelled_onset_chan', 'imprint_chan','EI_chan', 'PLHG_chan',...
        'resected_chan', 'labelled_onset_roi','resected_roi', ...
        'Surgery_outcome', 'Surgery_year', 'Outcome_year', 'Op_type'};
    onset_output = cell(size(sz_count_pat,1), 15);
    onset_output = cell2table(onset_output, 'VariableNames', colNames);
    % Populate patient Ids
    onset_output.Patient_id = string(sz_count_pat(:,1));
    
    % We also want to store WHEN onset is captured (to determine if there are any large disparities
    % For EI, I will take the mean if there is a variety of times
    
    %% Iterate through included patients, calculate onset and store in table
    % Window size for confirming seizure activity based on imprint
    wind_size = 2;
    
    incl_patients = onset_output.Patient_id;
    for pat = 1:size(onset_output,1)
        patient = incl_patients(pat);
        pat_seizures = string(data_tbl.patient_id) == patient;
        pat_data = data_tbl(string(data_tbl.patient_id) == patient,:);
        pat_imprints = cell_imprint(pat_seizures,:);
        sz_count = size(pat_imprints,1);
        
        % Calculate onset based on EI 
        [EI_tbl,Nd_tbl,~,~] = ms_epi_ind(pat_data);
        tbl_plhg = table();
        
        channel_details = json_data(1).channel_details;
        recorded_channels = channel_details(ismember(strrep(channel_details.chan_name, ' ', ''),...
            strrep(data_tbl.segment_channel_labels{1},' ', '')),:);

        recorded_roi_all_atlas = cat(2,recorded_channels.ROIname{:});
        incl_roi = recorded_roi_all_atlas(3,:)';
        unq_roi = unique(incl_roi, 'stable');

        n_chan =  size(pat_data.segment_data{1,1},1);
        
        % Create matrices to store onset channels for each method
        onset_imprint = zeros(n_chan, sz_count);
        onset_ei = onset_imprint;
        onset_plhg = onset_imprint;
        onset_time_imprint = zeros(1,sz_count);
       
        % Iterate over all patient seizures
        for sz = 1:sz_count
            % Onset based on seizure imprint
            sz_imprint = pat_imprints{sz,1}{1,1}; %.cell_imprint{1,1}
            imprint_wind = zeros(size(sz_imprint,1),size(sz_imprint,2)-(wind_size-1));
            for epoch = 1:(size(sz_imprint,2)-(wind_size-1))
                imprint_wind(:,epoch) = sum(sz_imprint(:,epoch:epoch+(wind_size-1)),2);
            end
            imprint_wind = imprint_wind==wind_size;
            chan_sum = sum(imprint_wind);
            % Identify when first activity is found
            onset_time = find(chan_sum,1,'first');
           
            if isempty(onset_time)
                sprintf("Patient %s, seizure %s: No imprint onset detected", patient, string(pat_imprints.segment_id(sz)))
                onset_time_imprint(sz) = NaN;
            else
                onset_imprint(:,sz) = imprint_wind(:,onset_time);
                onset_time_imprint(sz) = onset_time; %Imprint
            end
            
            % Onset based on Epileptigonicity Index
            %[EI,Nd,time,ERmaster,Na] = epileptogenicityIndex(ict');
            onset_ei(EI_tbl{sz}>0.3,sz) = 1;
        
            % Onset based on Phase-locked high-gamma
            [tbl_plhg(sz,:)] = ms_PLHG(pat_data(sz,:));
            [~,plhg_id] = mink(tbl_plhg(sz,:).when_plhg,4);
            onset_plhg(plhg_id,sz) = 1;
             
        end
       
        soz_roi = zeros(length(unq_roi),1);
        resected_roi = zeros(length(unq_roi),1);
        
        for grp = 1:length(unq_roi)
            soz_roi(grp) = sum(recorded_channels.is_soz(string(incl_roi) == string(unq_roi(grp))))>0;
            resected_roi(grp) = sum(recorded_channels.is_resected5(string(incl_roi) == string(unq_roi(grp))))>0;
        end
        
        % Store segment ids in output table (for future reference)
        onset_output(pat,:).Segment_ids = mat2cell(pat_data.segment_id',1,size(pat_data,1));
        % Store onset results in output table
        onset_output(pat,:).channel_names = {recorded_channels.chan_name};
        onset_output(pat,:).ROI_ids = mat2cell(unq_roi,length(unq_roi),1);
        onset_output(pat,:).labelled_onset_chan = {recorded_channels.is_soz};
        onset_output(pat,:).imprint_chan = mat2cell(onset_imprint,n_chan,sz_count);
        onset_output(pat,:).EI_chan = mat2cell(onset_ei,n_chan,sz_count);
        onset_output(pat,:).PLHG_chan = mat2cell(onset_plhg,n_chan,sz_count);
        onset_output(pat,:).resected_chan = {recorded_channels.is_resected5};
        onset_output(pat,:).labelled_onset_roi = mat2cell(soz_roi, length(unq_roi), 1);
        onset_output(pat,:).resected_roi = mat2cell(resected_roi,length(unq_roi),1);
        onset_output(pat,:).Surgery_outcome = pat_data.ilae(1);
        surgery_date = getfield(json_data(1).treatment_details.treatment_date, 'x_date');
        onset_output(pat,:).Surgery_year = cellstr(string(surgery_date(1:4)));
        onset_output(pat,:).Outcome_year = mat2cell(pat_data.ilae_year{1,1}, length(pat_data.ilae_year{1,1}), 1);
        onset_output(pat,:).Op_type = cellstr(pat_data.op_type(1));
    
        % Store when onset was detected for all methods
        when_onset_ei = zeros(1, sz_count);
        for sz = 1:sz_count
           when_onset_ei(sz) = Nd_tbl{sz}(EI_tbl{sz}==1);
        end
        when_onset = [onset_time_imprint' when_onset_ei'  min(tbl_plhg.when_plhg,[],2,'omitnan')];
        onset_output.when_onset{pat} = when_onset;
     

    end
    
end

