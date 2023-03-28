function [onset_output] = onset_chan_to_roi(data_tbl, json_data, onset_output, opts)
% Compute average EEG across regions of interest for each patient

% input:
%   - data_tbl: data table for ONE patient
%   - json_data: json data for same patient
%   - Onset_output: output table with onset detected on a channel level for
%   this patient
%   - optional inputs
%       - method: method of automatic detection used ("imprint", "EI", or
%       "PLHG")

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl
        json_data
        onset_output
        opts.method (1,1) string {mustBeMember(opts.method, ["imprint", "EI", "PLHG"])} = "imprint" % state method of automatic onset detection
    end
    
    %fill in optional arguments
    method = opts.method;

    if method == "imprint"
        onset = onset_output.Imprint_onset;
    elseif method == "EI"
        onset = onset_output.EI_onset;
    elseif method == "PLHG"
        onset = onset_output.PLHG_onset;
    end
  
    % Add empty column to store output
    onset_output = [onset_output cell(size(onset_output,1),1)];
    onset_output.Properties.VariableNames(size(onset_output,2)) = sprintf("%s_roi", method);
    
    pat_onset = onset{1,1};
    unq_roi = onset_output.ROI_ids{1,1};

    channel_details = json_data(1).channel_details;
    recorded_channels = channel_details(contains(strrep(channel_details.chan_name, ' ', ''),...
        strrep(data_tbl.segment_channel_labels{1},' ', '')),:);
    recorded_roi_all_atlas = cat(2,recorded_channels.ROIname{:});
    incl_roi = recorded_roi_all_atlas(3,:)';
    onset_roi = zeros(length(unq_roi), size(data_tbl,1));

    for grp = 1:length(unq_roi)
        onset_roi(grp,:) = sum(pat_onset(string(incl_roi) == string(unq_roi(grp)),:),1);
    end
    
    onset_roi(onset_roi>0) = 1;
    onset_output(:,size(onset_output,2)) =  mat2cell(onset_roi,size(onset_roi,1),size(onset_roi,2));

end