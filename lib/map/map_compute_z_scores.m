function roi_z = map_compute_z_scores(roi_data,norm_mean,norm_std)
% MAP_COMPUTE_Z_SCORES Compute z-scores of ROI data using normative map as
% a reference.
%
%   roi_z = MAP_COMPUTE_Z_SCORES(roi_data,norm_mean,norm_std) returns the
%   z-scored ROI data, roi_z, from segment ROI data (roi_data, ROIs x
%   measure features x iEEG segment) z-scored using the normative map mean
%   (norm_mean, ROIs x measure features) and standard deviation (norm_std,
%   ROIs x measure features).
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    roi_data (:,:,:) {mustBeNumeric}
    norm_mean (:,:) {mustBeNumeric}
    norm_std (:,:) {mustBeNumeric}
end
    
% check for standard deviation of zero
if sum(sum(norm_std==0)) > 0
    warning('Normative map contains values with standard deviation; will produce Inf z-scores')
end

% sizes
n_roi = size(roi_data,1);
n_feat = size(roi_data,2);
n_segm = size(roi_data,3);

% check size of normative map matches ROI data
if size(norm_mean,1) ~= n_roi || size(norm_std,1) ~= n_roi
    error('Number of ROIs in normative map does not match number in ROI data to be z-scored')
end
if size(norm_mean,2) ~= n_feat || size(norm_std,2) ~= n_feat
    error('Number of features in normative map does not match number in ROI data to be z-scored')
end

% tile normative map (explicitly tiles to avoid any possible confusion if
% dimension sizes are equal)
norm_mean_tile = repmat(norm_mean,1,1,n_segm);
norm_std_tile = repmat(norm_std,1,1,n_segm);

% z-score
roi_z = (roi_data - norm_mean_tile)./norm_std_tile;







