function [rel_bp,band_settings] = meas_rel_bp(pxx,freq,freq_bands,opts)
% MEAS_REL_BP Compute relative log10 bad power of each channel from tbe 
% power spectal densities.
%
% Also returns a structure of the band power settings.
% 
%   INPUTS
%
%       pxx: matrix of channel power spectal densities (channels x
%       frequencies)
%
%       freq: vector of frequencies in pxx
%
%       freq_bands: array of the lower (first column) and upper (second
%       column) bounds of the frequency bands to compute
%
%       Optional Name/Value pair arguments:
%           - "CombineBands" is a numeric vector specifying which bands in
%           freq_bands contribute to the final frequency bands. Each entry
%           in the vector corresponds to one row in freq_bands.
%           E.g., 
%           freq_bands = [1 4; 4 8; 8 10; 15 20]
%           combine_bands = [1 2 3 3]
%           rel_bp = meas_rel_bp(pxx,freq,freq_bands,"CombineBands",combine_bands)
%           returns a n channel x 3 array of relative log10 band power  in
%           1-4 Hz
%           4-8 Hz
%           8-10 + 15-20 Hz
%
%           The length of the CombineBands vector must be the same as the
%           number of rows in freq_bands, and it must not skip any
%           integers.
%
%           By default, no frequency bands are combined. 
%
% See also MEAS_PWELCH_PSDS, MEAS_ALL_PWELCH_PSD, MEAS_ALL_FROM_PSD
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    pxx {mustBeNumeric}
    freq {mustBeNumeric}
    freq_bands (:,2) {mustBeNumeric}
    opts.CombineBands = [];
end

% check pxx and freq dimensions)
if size(pxx,2) ~= length(freq)
    error('Number of columns in pxx must match the length of freq')
end

% number of channels
n_chan = size(pxx,1);

% number of frequency bands
n_bands = size(freq_bands,1);

% initialise array for storing band power
bp = zeros(n_chan,n_bands);

% compute band power from PSD
for b=1:n_bands
    bp(:,b) = bandpower(pxx',freq,freq_bands(b,:),'psd');
end

% combine frequency bands if requested
if ~isempty(opts.CombineBands)
    % new number of bands
    n_bands_new = max(opts.CombineBands);
    
    % check CombineBands format
    if ~isequal(1:n_bands_new,unique(opts.CombineBands))
        error('CombineBands vector must not skip integers')
    end
    if length(opts.CombineBands) ~= n_bands
        error('length of CombineBands must be the same as the number of rows in freq_bands')
    end
    
    % initialise array for combined band power
    bp_combined = zeros(n_chan,n_bands_new);
    
    % combine requested columns of bp
    for b=1:n_bands_new
        bp_combined(:,b) = sum(bp(:,opts.CombineBands==b),2);
    end
    
    % replace with combined band power
    bp = bp_combined;
    n_bands = n_bands_new;
    
else
    opts.CombineBands = 1:n_bands;
end

% save band settings
band_settings = struct();
band_settings.freq_bands = freq_bands;
band_settings.combine_bands = opts.CombineBands;

% log transform
bp = log10(1+bp);

% relative band power
rel_bp = bp./repmat(sum(bp,2),1,n_bands);


