function ieeg_data_reref = pre_reref(ieeg_data,ref)
% PRE_REFEF  Re-references iEEG data to a common average reference (mean or
% median average).
%
%   ieeg_data_reref = PRE_REREF(ieeg_data,ref) applies the specified
%   type of reference, ref, to channels x time iEEG data, ieeg_data. 
%
%   INPUTS:
%
%       ieeg_data: matrix of electrode signals (# channels x # time points)
%
%       ref: string specifying type of rereferencing; options are
%       (case sensitive)
%           -'CAR': common average reference (mean signal subtracted)
%           -'medianCAR': common average reference (median signal
%           substracted, useful if some channels have noise/outliers and
%           those channels cannot be removed)
%       For CAR and medianCAR referencing options, NaNs in a subset of 
%       channels are ignored when computing the average (i.e., the average
%       is computed from non-NaN channels).
%
%   OUTPUTS:
%
%       ieeg_data_reref: matrix of rereferenced iEEG data, same size as
%       ieeg_data
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022
% Modified from preproc_reref (GMS, 12 July 2021)
% Can be modified in future to add additional rereferencing options

% argument validation
arguments
    ieeg_data {mustBeNumeric}
    ref char {mustBeMember(ref,{'CAR','medianCAR'})}
end

% number of channels
n_chan=size(ieeg_data,1);

% rereference
switch ref
    
    % common average reference (mean)
    case 'CAR'
        average_signal = mean(ieeg_data,1,'omitnan');                       % mean signal
        ieeg_data_reref = ieeg_data - repmat(average_signal,[n_chan, 1]);   % remove mean signal
    
    % common average reference (median)
    case 'medianCAR'
        average_signal = median(ieeg_data,1,'omitnan');                     % median signal
        ieeg_data_reref = ieeg_data - repmat(average_signal,[n_chan, 1]);   % remove median signal    
        
end

      

