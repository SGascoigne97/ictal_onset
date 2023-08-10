function f = vis_bad_chan(json_data,segm_idx,mat_path,options)
% VIS_BAD_CHAN Plots the iEEG traces and PSDs of one iEEG segment, with bad
% channels from each category plotted in a different colour than "good"
% channels. 
%
% Bad channels must be identified first and added to the JSON structure
% (see functions in "see also" section).
%
%   f = VIS_BAD_CHAN(json_data,segm_idx,mat_path) returns the figure handle
%   f to the bad channel visualisation of segment number segm_idx in the
%   JSON structure json_data. The iEEG workspaces for the corresponding
%   database export must be in the directory mat_path. 
%
%   See argument validation notes for optional name/value pair arguments. 
%   In particular, the bad channel labels plotted can be modified using
%   'BadChanFields'.
%
%   See also VIS_EEG, VIS_PSD, PRE_FILT_CHAN_BY_FEATURE, PRE_DETECT_BAD_CHAN, 
%   PRE_ALL_DETECT_BAD_CHAN, PRE_ADD_BAD_CHAN_TO_JSON, PRE_FIND_MISS_CHAN,
%   PRE_FIND_CHAN_WITH_NO_MAPPING, PRE_CONSOLIDATE_BAD_CHAN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    json_data struct
    segm_idx (1,1) {mustBeInteger}
    mat_path char

    % Name/Value pairs

    % Bad channel fields to use for bad channel names; will plot one
    % subplot for each field
    options.BadChanFields cell = {'bad_chan_filtered',...
        'bad_chan_detected','bad_chan_mat','bad_chan_miss','bad_chan_no_mapping',...
        'bad_chan_all'};

    % Colors
    options.BadColor (1,3) {mustBeNumeric} = [246,80,80]/255; % bad channels color
    options.NotBadColor (1,3) {mustBeNumeric} = [1 1 1]*0.75; % not bad channels color
    options.TintOrShade (1,1)...                                    % for every other channel colour, lightens (tints) if < 0 and darkens (shades) if > 0
        {mustBeNumeric, mustBeLessThan(options.TintOrShade,1),...   % transformed so closer to 1/-1 = more impact
        mustBeGreaterThan(options.TintOrShade,-1)} = 0.2;           % only used for iEEG traces
    
    % Figure/plot settings
    options.FigPos (1,4) {mustBeNumeric} = [2 2 50 30];       % figure size
    options.Offset (1,1) {mustBeNumeric} = 200;               % spacing between traces

    % PSD settings
    options.Window (1,1) {mustBeNumeric} = 2;         % pwelch window size in seconds (see meas_pwelch_psd)
    options.Overlap (1,1) {mustBeNumeric} = 1;        % pwelch window overlap size in seconds (see meas_pwelch_psd)
    options.Freq (1,:) {mustBeNumeric} = 0.5:0.5:100; % frequencies at which to compute PSD

end

% check that all requested bad channel fields are in the JSON structure
json_fields = fieldnames(json_data);                                % json structure field names
miss_fields = setdiff(options.BadChanFields,json_fields,'stable');  % any fields not in json structure?
if ~isempty(miss_fields)
    miss_fields = cellfun(@string,miss_fields) + ' ';
    warning(['Following fields are missing and will not be plotted: ' miss_fields{:}])
end
options.BadChanFields = intersect(options.BadChanFields,json_fields,'stable'); % only use fields present in both

% number of bad channel fields
n_fields = length(options.BadChanFields);

% load segment
eeg_fn = db_all_json_eeg_fn(json_data);
eeg_fn = eeg_fn{segm_idx};
load([mat_path '/' eeg_fn],'eeg_channels','eeg_data','eeg_fs')
n_chan = length(eeg_channels);

% booleans of which channels are bad in each category
bad_chan_bool = zeros(n_chan,n_fields);
for i=1:n_fields

    % get names of bad channels
    bad_chan = json_data(segm_idx).(options.BadChanFields{i});

    % find bad channels in segment 
    bad_chan_bool(:,i) = ismember(eeg_channels,bad_chan);
    
end
bad_chan_bool = logical(bad_chan_bool);

% plot
f=figure();
set(f,'units','centimeters','position',options.FigPos)
n_row = 3;
for i=1:n_fields

    % bad channels according to field
    is_bad = bad_chan_bool(:,i);

    % color array 
    clrs = repmat(options.NotBadColor,n_chan,1);
    clrs(is_bad,:) = repmat(options.BadColor,sum(is_bad),1);

    % plot eeg
    if i==1
        ylab = 'channel';
    else
        ylab = '';
    end
    subplot(n_row,n_fields,[i i+n_fields])
    vis_eeg(eeg_data,eeg_fs,'PlotNewFig',false,'Color',clrs,...
        'Offset',options.Offset,'ChannelNames',eeg_channels,...
        'TintOrShade',options.TintOrShade,'YAxisLabel',ylab);
    title(options.BadChanFields{i},'Interpreter','none')

    % plot psd
    if i==1
        ylab = 'PSD (dB/Hz)';
    else
        ylab = '';
    end
    subplot(n_row,n_fields,i+n_fields*2)
    vis_psd(eeg_data,eeg_fs,'EEG',...
        'Window',options.Window,'Overlap',options.Overlap,'Freq',options.Freq,...
        'Color',clrs,'PlotNewFig',false,'YAxisLabel',ylab);

end

% title with segment info
patID = json_data(segm_idx).patient_details.patID;
hosp = json_data(segm_idx).patient_details.Hospital;
segm_num = db_all_segm_num(json_data);
segm_num = segm_num(segm_idx);
sgtitle([hosp ' ' patID ', segment ' num2str(segm_num)])



