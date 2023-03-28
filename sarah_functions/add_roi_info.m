function [data_tbl_new] = add_roi_info(data_tbl, opts)
% Add region of interest (ROI) information for patients in data table
% 
% input:
%   - data_tbl: full data table

% output
%   - data_tbl_new: full data table with additional ROI information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl
        opts.roi_parc (1,1) double {mustBeNumeric} = 3 % Chosen parcellation scheme for ROIs (between 1 and 4)
    end
    
    %fill in optional arguments
    roi_parc = opts.roi_parc;

    % Extract patient IDs
    patients = unique(data_tbl.patient_id);

    %% We will use ROIs to identify the seizure onset zone
    % Add column to save ROI information
    data_tbl_new = data_tbl;
    data_tbl_new.channel_roi = cell(size(data_tbl_new,1),1);
    % Import channel information
    
    channel_folder = '../../channels_ROI/';
    % Add information to data table (ROI, is channel identified as SOZ, and
    % were channels resected at 3 and 5 mm)
    for pat=1:(length(patients))
        % Extract seizure data for this patient
        pat_data = data_tbl_new(string(data_tbl_new.patient_id) == patients{pat},:);
        % Import channel information (including ROI identification) for this
        % patient
        channel_info = load([channel_folder patients{pat} '/channels.mat']);
        % Identify channels included in EEG recordigns (some channels were
        % removed during data preprocessing
        incl_chan = ismember(string(channel_info.channels.name), string(pat_data.segment_channel_labels{1,1}));
         try
             % Extract ROI names for each channel
             ROI_incl_channels = channel_info.channels.ROIname(incl_chan,roi_parc);
             % Add ROI information to data table
            data_tbl_new.channel_roi(string(data_tbl_new.patient_id) == patients{pat} ) = {ROI_incl_channels};
         catch
            warning('Patient %s: ROI names not available', patients{pat});
        end
       
        
        % Also add a binary array identifying which channels are deemed to be
        % seizure onset by clinicians
        data_tbl_new.channel_soz(string(data_tbl_new.patient_id) == patients{pat} ) = {channel_info.channels.is_soz(incl_chan,:)};
        try
            % Add resection information for each channel - this is not
            % available foir all patients. A warning will be printed for
            % patints lacking resection information
            data_tbl_new.channel_resected_3mm(string(data_tbl_new.patient_id) == patients{pat} ) = {channel_info.channels.isResected3mm(incl_chan,:)};
            data_tbl_new.channel_resected_5mm(string(data_tbl_new.patient_id) == patients{pat} ) = {channel_info.channels.isResected5mm(incl_chan,:)};
        catch
            warning('Patient %s: Resection info not available', patients{pat});
        end
    end
end