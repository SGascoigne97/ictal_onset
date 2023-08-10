function f = vis_noise_detec_manual(pat_json, preproc_label, path_pipeline_exports, opts)
% VIS_NOISE_DETEC_MANUAL Creates plots for manual detection of noise in the
% preictal segment
%
% Sarah J. Gascoigne
% CNNP Lab, Newcastle University
% March 2023

 arguments
        pat_json struct
        preproc_label string
        path_pipeline_exports string
        % optional
        opts.offset double = 400
 end
    offset = opts.offset;
    % Add detais on how the function works
    
    % Extract patient ID
    pat = string(db_extract_json_metadata(json_data(1),'patID'));
    % Obtain file path for all patient seizures 
    pat_json_preproc_fn = db_all_json_preproc_fn(pat_json, preproc_label);
    % Extract the unique bad channels for the patient
    bad_chan = [];
    for sz = 1:size(pat_json,1)
        bad_chan = [bad_chan; string(pat_json(sz).bad_chan_all)] ;
    end
    bad_chan_unq = unique(bad_chan);
    
    %% Select a seizure and extract the EEG data and PSD
    for sz = 1:size(pat_json,1)
        
        raw_data_location = string(pat_json(sz).eeg_fn);
        raw_data_struct = load(raw_data_location);
        raw_data = raw_data_struct.eeg_data;
        raw_channels = raw_data_struct.eeg_channels;
    
        preproc_data_location = string(['export_abnorm/', pat_json(sz).fn_analysis.preEurope.preproc]);
        preproc_data_struct = load(preproc_data_location);
        preproc_data = preproc_data_struct.eeg_data;
        preproc_channels = preproc_data_struct.eeg_channels;

        % load EEG data
        load([path_pipeline_exports '/' pat_json_preproc_fn{sz}],...
            'eeg_data','eeg_channels','eeg_fs','x_oid');
        
        % Plot EEG data for selected seizure
        % Set colours to grey except for removed channels
        colours = ones(length(raw_channels), 3)*0.7;
        % Add highlights for removed channels
        colours(contains(raw_channels, bad_chan_unq),:) = repmat([1,0,0],...
            sum(contains(raw_channels, bad_chan_unq)), 1);
        f = figure((sz*2)-1);
        set(f,'units','centimeters','position',[2 2 20 20]);
        subplot(121)
        vis_eeg(raw_data, eeg_fs, 'Color', colours, 'Offset', offset,...
            'ChannelNames',raw_channels, 'PlotNewFig', false);
        xlim([0, pat_json(sz).eeg_pre])
        xline(110, 'k', 'LineWidth', 2, 'LineStyle','--')
        title('Raw EEG with noisy channels highlighted')
    
        subplot(122)
        vis_eeg(preproc_data, eeg_fs, 'Offset', offset,...
            'ChannelNames',preproc_channels, 'PlotNewFig', false);
        xlim([0, pat_json(sz).eeg_pre])
        xline(110, 'k', 'LineWidth', 2, 'LineStyle','--')
        title('Preprocessed EEG')
        sgtitle(sprintf('Patient %s, seizure %s \n EEGs', pat, x_oid))
    
        %Plot PSD for preictal segment (for visual idnetification of noisy
        %channels)
    
        % extract pre-ictal segment (excl 10 seconds closest to seizure onset)
        pre_data = eeg_data(:,1:(110*eeg_fs));
        vis_psd(pre_data,eeg_fs,'EEG', 'Interactive', true,...
            'ChannelNames',eeg_channels);
        title(sprintf('Patient %s, seizure %s \n PSD of preprocessed EEG', pat, x_oid))
        
        
    end
end