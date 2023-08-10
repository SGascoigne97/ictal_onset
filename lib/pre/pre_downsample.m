function ieeg_data_proc = pre_downsample(ieeg_data,fs,fs_final,opts)
% PRE_DOWNSAMPLE Downsample iEEG segment to the desired frequency.
%
%   ieeg_data_proc = PRE_DOWNSAMPLE(ieeg_data,fs,fs_final) downsamples the
%   iEEG segment ieeg_data (channels x time numeric array) from sampling 
%   frequency fs Hz to sampling frequency fs_final Hz. The downsampled iEEG
%   segment is returned as ieeg_data_proc (channels x time numeric array).
%
%   ieeg_data_proc = PRE_DOWNSAMPLE(...,Name, Value) specifies the
%   following options for the segment downsampling:
%       - 'Method' (string) specifies the method to use for downsampling.
%       Current options are 'Interp' (default), which uses linear
%       interpolation (MATLAB function interp1)
%
% See also INTERP1, RESAMPLE
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022
%
% TODO: add "Resample" as 'Method' option


arguments
    ieeg_data {mustBeNumeric}
    fs (1,1) {mustBeNumeric}
    fs_final (1,1) {mustBeNumeric}
    opts.Method char {mustBeMember(opts.Method,{'Interp'})} = 'Interp'
end

% check that fs_final is lower than fs
if fs_final > fs
    error('Final sampling frequency must be lower than original sampling frequency')
end

switch opts.Method
    % linear interpolation
    case 'Interp'
        % current time vector
        n = size(ieeg_data,2);
        t = (1:n)/fs;
        
        % new time vector
        t_final = (1:floor(t(end)*fs_final))/fs_final;
        
        % linear interpolation
        ieeg_data_proc = interp1(t,ieeg_data',t_final);
        ieeg_data_proc = ieeg_data_proc';
end