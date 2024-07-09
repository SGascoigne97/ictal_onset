function [f,h] = vis_psd(data,freq_info,data_type,options)
% VIS_PSD Plot PSDs from iEEG data or pre-computed PSDs. Can also make
% interactive (bolds PSD from one channel at a time when clicked).
%
%
%   [f,h] = VIS_PSD(eeg_data,fs,'EEG') computes and plots PSDs from the
%   iEEG data eeg_data (channels x time) with sampling frequency fs Hz. The
%   figure handle f and plot handle h are returned.
%
%   [f,h] = VIS_PSD(psd,freq,'PSD') plots the PSDs in psd (channels x 
%   frequency) with the corresponding frequency vector freq.
%
%   [f,h] = VIS_PSD([],[],'EEG','Filename',eeg_fn) computes and plots 
%   PSDs from the iEEG data stored in the mat file located at eeg_fn (must  
%   be full path to file from current location). The mat file must contain  
%   the variables eeg_data, eeg_fs, and eeg_channels. The channel PSDs will  
%   be automatically labelled with a legend using this syntax.
%
%   [f,h] = VIS_PSD([],[],'PSD','Filename',psd_fn) computes and plots the
%   PSDs stored in the mat file located at psd_fn (must be full path to
%   file from current location). The mat file must contain the variables
%   pxx, freq, and eeg_channels. The channel PSDs will be automatically
%   labelled with a legend using this syntax.
%
%
%   See argument validation for optional name/value pair arguments. Use
%   'Interactive' to make plot interactive.
%
%   See also MEAS_PWELCH_PSD, VIS_EEG, VIS_PSD_SUMMARIES
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023; modified from earlier version

arguments
    data {mustBeNumeric}
    freq_info {mustBeNumeric}
    data_type char {mustBeMember(data_type,{'EEG','PSD'})}
    
    % Filename for EEG/PSD data if data and freq_info are empty
    options.Filename char = '';
    
    % PSD settings - only used if data_type = 'EEG'
    options.Window = 2;         % pwelch window size in seconds (see meas_pwelch_psd)
    options.Overlap = 1;        % pwelch window overlap size in seconds (see meas_pwelch_psd)
    options.Freq = 0.5:0.5:100; % frequencies at which to compute PSD
    
    % Controls colours
    % Colormap/ColormapN are ignored if Color is provided
    options.Color {mustBeNumeric} = [];             % colour(s) to use for lines - if < # of channels, will repeat
    options.Colormap char = 'lines';                % if options.Color is empty, name of colourmap to use
    options.ColormapN (1,1) = 7;                    % if colours set using colormap, number of colours to draw (will repeat to match channel number)
    
    % Controls labels
    options.ChannelNames cell = {};                     % cell array of channel names to use to label PSDs; used to plot legend if not empty
    options.ForceLegendOff = false;                     % if true, turns off channel legend even if ChannelNames is not empty
    options.XAxisLabel char = 'frequency (Hz)';         % x-axis label
    options.YAxisLabel char = 'PSD (dB/Hz)';            % y-axis label
    options.Title char = '';                            % title
    options.FontSize (1,1) {mustBeNumeric} = 12;        % font size for tick labels
    options.LabelFontSize (1,1) {mustBeNumeric} = 14;   % font size for axis labels and title
    options.NLegendCol (1,1) {mustBeNumeric} = 2;       % number of columns in legend
    options.LegendFontSize (1,1) {mustBeNumeric} = 8;   % font size for legend
    options.LegendLocation char = 'westoutside';        % legend location
    % Controls other figure/plot settings
    options.Interactive {mustBeNumericOrLogical} = false;   % whether to make plot interactive (can click on and bold PSD line)
    options.PlotNewFig {mustBeNumericOrLogical} = true;     % whether plot new figure (and pass figure handle f as output); turn off to plot as subplot in existing figure
    options.FigPos = [2 2 30 20];                   % figure position in centimeters; only used for new figure
    options.Linewidth (1,1) = 1;                    % linewidth of lines used for iEEG traces
    options.AxisLim = [];                           % axis limits - should be (1,4) or axis style (e.g., 'square')
end

%%% COMPUTE/GET PSD
% load from filename if variables are empty
if isempty(data) && isempty(freq_info)
    
    switch data_type
        
        % EEG case
        case 'EEG'
            load(options.Filename,'eeg_data','eeg_fs','eeg_channels'); % load
            data = eeg_data; freq_info = eeg_fs; clearvars eeg_data eeg_fs
            
            % compute
            [pxx,freq] = meas_pwelch_psd(data,options.Window,options.Overlap,options.Freq,freq_info);
            
        % PSD case
        case 'PSD'
            load(options.Filename,'pxx','freq','eeg_channels'); % load
    end
    
    % channels
    if ~isempty(options.ChannelNames)
        warning('Provided channel name will be overwritten by eeg_channels variable in mat file')
    end
    options.ChannelNames = eeg_channels; clearvars eeg_channels
    
% filename should not be provided if eeg_data and eeg_fs are not both empty
elseif ~isempty(options.Filename)
    error('data and freq_info should be empty if Filename is provided')

% if PSD data provided directly, change names 
elseif isequal(data_type,'PSD')
    
    % check sizes
    if ~isequal(size(data,2),length(freq_info))
        error('The number of columns in data should be equal to the length of freq_info when data_type is PSD')
    end
    
    % rename
    pxx = data;
    freq = freq_info;
    
% if EEG data is provided directly, compute PSD variables
elseif isequal(data_type,'EEG')
    
    % check size of freq_info (should be scalar)
    if length(freq_info) ~= 1
        error('freq_info should be a scalar of the EEG sampling frequency if data_type is EEG')
    end
    
    % compute
    [pxx,freq] = meas_pwelch_psd(data,options.Window,options.Overlap,options.Freq,freq_info);

end
clearvars data freq_info % remove old non-specific variables

% number of channels
n_chan = size(pxx,1);


%%% SET COLOURS

% if Color is empty, set using Colormap name + requested number of colors 
% to draw from heatmap (ColormapN)
if isempty(options.Color)
    eval(['options.Color = ' options.Colormap '(' num2str(options.ColormapN) ');']);
% otherwise, check that Color has three columns
elseif size(options.Color,2) ~= 3
    error('Color must either be empty or an n x 3 array')
end

% if number of rows in Color is less than number of EEG channels, repeat
% rows of Color to match
options.Color = rep_colors(options.Color,n_chan);


%%% PLOT
% start new figure
if options.PlotNewFig
    f = figure();
    set(f,'units','centimeters','position',options.FigPos);
else
    f = gcf;
end
    
% plot
UD.h = plot(freq,10*log10(pxx'),'LineWidth',options.Linewidth);
set(gcf,'Color',[1 1 1]); % background colour
h=UD.h; % figure handle to return

% set colors
for i=1:n_chan
    set(UD.h(i),'Color',options.Color(i,:));
end

% axis limits
if ~isempty(options.AxisLim)
    axis(options.AxisLim)
else
    xlim([0 Inf])
end

% set font sizes and labels
set(gca,'FontSize',options.FontSize);
xlabel(options.XAxisLabel,'FontSize',options.LabelFontSize);
ylabel(options.YAxisLabel,'FontSize',options.LabelFontSize);
if ~isempty(options.Title)
    title(options.Title,'FontSize',options.LabelFontSize);
end

% legend
if ~isempty(options.ChannelNames) && ~options.ForceLegendOff
    legend(UD.h,options.ChannelNames,...
        'location',options.LegendLocation,'NumColumns',options.NLegendCol,...
        'FontSize',options.LegendFontSize);
end


%%% MAKE INTERACTIVE
if options.Interactive

    % can only have one interactive plot per figure, so enforce PlotNewFig = true
    if ~options.PlotNewFig
        error('Can only plot interactive plot if PlotNewFig = true')
    end
    
    % change color and line width when click
    set(UD.h,'ButtonDownFcn',{@line_click,f});
    UD.Selected = [];
    clrs_bold = get(UD.h,'Color');
    clrs_light = cell(length(clrs_bold),1);
    for i=1:length(clrs_bold)
        clrs_light{i} = clrs_bold{i} + ((1 - clrs_bold{i})*0.5);
        set(UD.h(i),'Color',clrs_light{i});
    end
    UD.clrs_light = clrs_light;
    UD.clrs_bold = clrs_bold;
    f.UserData = UD;

end
end

% clicking function
function line_click(obj_h,event_data,f)

UD = f.UserData;

n_lines = length(UD.h);
lw = zeros(n_lines,1);
for i=1:n_lines
    lw(i) = UD.h(i).LineWidth;
end
lw_light = min(lw);
lw_bold = lw_light*5;

clrs_light = UD.clrs_light;
clrs_bold = UD.clrs_bold;

idx = find(obj_h == UD.h);

if isequal(idx,UD.Selected)
    UD.Selected = [];
    UD.h(idx).LineWidth = lw_light;
    UD.h(idx).Color = clrs_light{idx};
else
    if ~isempty(UD.Selected)
        set(UD.h(UD.Selected),'Color',clrs_light{UD.Selected}); % reset color of previously selected
    end
    UD.Selected = idx; % new selected
    set(UD.h,'LineWidth',lw_light)
    set(UD.h(idx),'LineWidth',lw_bold);
    set(UD.h(idx),'Color',clrs_bold{idx});
    uistack(UD.h(idx),'top')
end

f.UserData = UD;
   
end

% function for repeating colors if fewer colours than objects to colour
% will also ensure final list of colors is same size (n rows) as number of objects 
function my_colors = rep_colors(my_colors,n_final)

n_row = size(my_colors,1);
if n_row < n_final
    my_colors = repmat(my_colors,ceil(n_final/n_row),1);
end
my_colors = my_colors(1:n_final,:); 

end