function [target_roi, target_roi_str] = chan_to_roi_crit(target_chan, incl_roi, unq_roi, opts)
% Function to convert channels to regions using inclusion criteria
% This can be used for converting onset or resection from channels to
% regions
%
%   >=Threshold% of channels are onset: region is included
%   0% of channels are onset: region is not included
%   <Threshold% of channels are onset: region is ignored

% Sarah J Gascoigne 24/04/2023

% input:
%   - target_chan: binary array idicating if channels are in target (e.g., onset/resected) 
%   - incl_roi: % cell array of the region label for each included channel
%   - optional inputs
%   	- threshold: proportion of channels in region to be included 


% output
%   - target_roi: binary array indicating if regions are in target (e.g., onset/resected) 

    arguments
        target_chan {mustBeNumericOrLogical}
        incl_roi % Cell array of the region label for each included channel
        unq_roi % Cell array of unique regions included
        opts.threshold (1,1) double = NaN % define as value within 0 and 1 (i.e. 0.25 = 25% of channels) 
                                          % Nan means that one channels in region satisfies criteria for inclusion
    end
    
    %fill in optional arguments
    threshold = opts.threshold;

    % Create empty array to store binary output (is region onset/resected)
    target_roi = zeros(1,length(unq_roi));

    for grp = 1:length(unq_roi)
        % Determine the number of channels in the region
        tot_chan = sum(ismember(incl_roi,unq_roi(grp)));
        % Determine the number of channels in the region which are onset
        total_target_chan = sum(target_chan(string(incl_roi) == string(unq_roi(grp))));
        
        % Implement inclusion criteria
        %   >=25% of channels are onset: region is onset
        %   0% of channels are onset: region is not onset
        %   <25% of channels are onset: region is ignored

        if isnan(threshold)
            min_chan = 1;
        else
            min_chan = ceil(tot_chan*threshold); % round up minimum channel count
        end
        
        % Store logical array indicating which regions are onset (ignored regions
        % set as NaN)
        if total_target_chan >= min_chan
            target_roi(grp) = 1;
        elseif total_target_chan == 0
            target_roi(grp) = 0;
        else
            target_roi(grp) = NaN;
        end
    end

    % Also store the name of regions (to check correct regions are
    % identified)
    unq_roi_str = string(unq_roi);
    target_roi(isnan(target_roi)) = 0;
    target_roi_str = unq_roi_str(find(target_roi));
