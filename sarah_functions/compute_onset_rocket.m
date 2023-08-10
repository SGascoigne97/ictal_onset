function [onset_imprint, onset_time_imprint] = compute_onset_rocket(pat_data, pat_imprints, opts)

% input:
%   - pat_data: % table with all seizures for subject
%   - pat_imprints: table of imprints for all seizures included in data table
%   - Optional arguments
%       - opts.wdw_size: Integer denoting window size for considering activity 
                       % as included in the onset (default 2 to reduce risk 
                       % of inclusion of spurious activity)
%       - opts.det: Integer denoting the number of seconds following onset 
                        % detection where activity will be included

% output
%   - onset_output: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
% pat_data

% Sarah Jane Gascoigne & Yujiang Wang
% 25/05/2023

    arguments
        pat_data 
        pat_imprints 
        opts.wdw_sz (1,1) double = 2 
        opts.det (1,1) double {mustBeNumeric} = 0
    end

    % Set optional arguments
    wdw_sz = opts.wdw_sz;
    det = opts.det;
    
    %% Iterate through seizures and compute onset for each method
    sz_count = size(pat_imprints,1); % Number of recorded seizures
    n_chan =  size(pat_data.segment_data{1,1},1); % Number of channels recorded

    % Create matrices to store onset channels for each method
    onset_imprint = zeros(n_chan, sz_count);
    onset_time_imprint = zeros(1,sz_count);
       
    % Iterate over all patient seizures
    for sz = 1:sz_count
        % Onset based on seizure imprint
        sz_imprint = pat_imprints{sz,1}{1,1}; 
        imprint_wind = zeros(size(sz_imprint,1),size(sz_imprint,2)-(wdw_sz-1));
        for epoch = 1:(size(sz_imprint,2)-(wdw_sz-1))
            imprint_wind(:,epoch) = sum(sz_imprint(:,epoch:epoch+(wdw_sz-1)),2);
        end
        imprint_wind = imprint_wind==wdw_sz;
        chan_sum = sum(imprint_wind);
        % Identify when first activity is found
        onset_time = find(chan_sum,1,'first');
       
        if isempty(onset_time)
            fprintf("Patient %s, seizure %s: No imprint onset detected \n", pat_data.patient_id{1}, string(pat_imprints.segment_id(sz)))
            onset_time_imprint(sz) = NaN;
        else
            onset_imprint(:,sz) = sum(imprint_wind(:,(onset_time:onset_time+det)),2) > 0;
            onset_time_imprint(sz) = onset_time; %Imprint
        end
        
       

    end
end