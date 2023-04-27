function [target_roi] = chan_to_roi_crit(target_chan, incl_roi, unq_roi, opts)
% Function to convert channels to regions using inclusion criteria
% This can be used for converting onset or resection from channels to
% regions
%
%   >=Threshold% of channels are onset: region is included
%   0% of channels are onset: region is not included
%   <Threshold% of channels are onset: region is ignored

% Sarah J Gascoigne 24/04/2023

% input:
%   - pat_onset: onset_output table for one specific patient 
%   - atlas_scale: table containing ROI atlas at chosen scale
%   - optional inputs
%       - det_method: state the method of onset detection you'd like to use
%           ("imprint", "EI", "PLHG")
%       - comparison: state whether the comparisons are based on
%           comparisons between onset and resected regions ("resction") or
%           between pairs of seizure onsets ("pairwise")

% output
%   - comp_table: table containing raw and normalised Hausdorff's distance

    arguments
        target_chan {mustBeNumericOrLogical}
        incl_roi % Cell array of the region label for each included channel
        unq_roi % Cell array of unique regions recorded 
        opts.threshold (1,1) double {mustBeInRange(opts.threshold, 0, 1)} = 0.25 % define as value within 0 and 1 (i.e. 0.25 = 25% of channels)
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
        min_chan = ceil(tot_chan*threshold); % round up minimum channel count
        
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
