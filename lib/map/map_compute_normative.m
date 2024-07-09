function [map_avg,map_spread,map_n] = map_compute_normative(roi_data_avg,opts)
% MAP_COMPUTE_NORMATIVE Compute normative map from ROI-level measure of
% multiple iEEG segments.
%
% INPUTS
%       
%   roi_data_avg: array of ROI-level data (averaged across channels), ROI x
%   measure features x iEEG segments. See outputs of map_all_chan2rois.
%
%   Optional arguments:
%       - 'AverageType': character specifying the type of centre and spread
%       measures; 'mean' = mean and standard deviation (default), 'median'
%       = median and median absolue deviation.
%       - 'MinN': minimum sample size (number of iEEG segments) for
%       normative map; data for any ROIs with smaller sample sizes will be
%       changed to NaNs. Default = 2 to avoid standard deviations of 0. 
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    roi_data_avg {mustBeNumeric}
    opts.AverageType char {mustBeMember(opts.AverageType,{'mean','median'})} = 'mean';  
    opts.MinN (1,1) = 2;
end

% sample size
not_nan = ~isnan(roi_data_avg);
map_n = sum(not_nan,3);

% dimension corresponding to iEEG segment
dim = 3;

% measure of average/centre and spread
switch opts.AverageType
    case 'mean'
        map_avg = mean(roi_data_avg,dim,'omitnan');         % mean
        map_spread = std(roi_data_avg,[],dim,'omitnan');    % standard deviation
    case 'median'
        map_avg = median(roi_data_avg,dim,'omitnan');       % median
        map_spread = mad(roi_data_avg,1,dim,'omitnan');     % median absolute deviation
end

% change to NaNs if ROI sample size less than min n
map_avg(map_n<opts.MinN) = NaN;
map_spread(map_n<opts.MinN) = NaN;
        
        
