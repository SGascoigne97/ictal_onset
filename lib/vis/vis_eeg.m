function [f,h] = vis_eeg(eeg_data,fs,options)
% VIS_EEG Plot the iEEG traces of one iEEG segment.
%
%   [f, h] = VIS_EEG(eeg_data,fs) plots the iEEG time series eeg_data
%   (channels x time) with sampling frequency fs Hz and returns the figure
%   and plot handles f and h. 
%
%   [f, h] = VIS_EEG([],[],'Filename',eeg_fn) plots the iEEG times series
%   stored in the mat file located at eeg_fn (must be full path to file
%   from current location). The mat file must contain the variables
%   eeg_data, eeg_fs, and eeg_channels. The channel time series will be 
%   labelled using this syntax.
%
%   See comments in argument validation for Name/Value pair arguments.
%
%   Has numerous optional arguments to control 
%   - spacing between traces ("Offset")
%   - labels ("TimeStart", "ChannelNames", "XAxisLabel", "YAxisLabel", 
%     "Title", "YFontSize", "XFontSize", "LabelFontSize", "ForceLabelsOff")
%   - colors ("Color", "Colormap", "ColormapN", "TintOrShade") 
%   - highlighted parts of iEEG ("HighlightTimes", "HighlightAlpha",
%   "HighlightColors")
%   - other aspects of the figure/plot ("PlotNewFig", "FigPos",
%   "Linewidth", "AxisLim")
%
%   OUTPUTS:
%
%       f: handle to figure if new figure requested (default); otherwise,
%       empty
%
%       h: handle to plot
%       
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023
%
% Modified from plotEEG (YW) and vis_plot_eeg (GMS/SG/YW)

arguments
    eeg_data {mustBeNumeric}
    fs {mustBeNumeric}
    
    % Filename for EEG data if eeg_data and fs are empty
    options.Filename char = '';
    
    % Controls spacing between traces
    options.Offset (1,1) {mustBeNumeric} = 200;     % size of gap between channel traces
    
    % Controls colours
    % Colormap/ColormapN are ignored if Color is provided
    options.Color {mustBeNumeric} = [];             % colour(s) to use for lines - if < # of channels, will repeat
    options.Colormap char = 'lines';                % if options.Color is empty, name of colourmap to use
    options.ColormapN (1,1) = 7;                    % if colours set using colormap, number of colours to draw (will repeat to match channel number)
    options.TintOrShade (1,1)...                    % for every other channel colour, lightens (tints) if < 0 and darkens (shades) if > 0
        {mustBeNumeric, mustBeLessThan(options.TintOrShade,1),... % transformed so closer to 1/-1 = more impact
        mustBeGreaterThan(options.TintOrShade,-1)} = 0; 
    
    % Controls highlight
    options.HighlightTimes (:,2) = [];              % beginning (1st column) and end (2nd column) times of each highlight; number of rows = number of highlights
    options.HighlightAlpha (1,1) = 0.5;             % transparancy of highlight
    options.HighlightColors (:,3) = [0.5 0.5 0.5];  % highlight colors; will repeat if < number of highlights
    
    % Controls labels
    options.TimeStart {mustBeNumeric} = [];         % where to start time vector for x-axis label; if is empty, will default to 1/fs later in code
    options.ChannelNames cell = {};                 % cell array of channel names to use to label y-axis
    options.ForceLabelsOff = false;                 % if true, turns off y-axis channel labels even if ChannelNames isn't empty
    options.XAxisLabel char = 'time (s)';           % x-axis label
    options.YAxisLabel char = 'channel';            % y-axis label
    options.Title char = '';                        % title
    options.YFontSize (1,1) {mustBeNumeric} = 6;        % font size for y-axis tick labels
    options.XFontSize (1,1) {mustBeNumeric} = 12;       % font size for x-axis tick labels
    options.LabelFontSize (1,1) {mustBeNumeric} = 14;   % font size for axis labels and title
    
    % Controls other figure/plot settings
    options.PlotNewFig {mustBeNumericOrLogical} = true; % whether plot new figure (and pass figure handle f as output); turn off to plot as subplot in existing figure
    options.FigPos = [2 2 20 20];                   % figure position in centimeters; only used for new figure
    options.Linewidth (1,1) = 1;                    % linewidth of lines used for iEEG traces
    options.AxisLim = [];                           % axis limits - should be (1,4) or axis style (e.g., 'square')
end


% load from filename if eeg_data and eeg_fs are empty
if isempty(eeg_data) && isempty(fs)
    
    % load data
    load(options.Filename,'eeg_data','eeg_fs','eeg_channels')
    fs = eeg_fs; clearvars eeg_fs
    
    % channels
    if ~isempty(options.ChannelNames)
        warning('Provided channel name will be overwritten by eeg_channels variable in mat file')
    end
    options.ChannelNames = eeg_channels; clearvars eeg_channels
    
% filename should not be provided if eeg_data and eeg_fs are not both empty
elseif ~isempty(options.Filename)
    error('eeg_data and fs should be empty if Filename is provided')
end

% set first time point to 1/fs if not specified
if isempty(options.TimeStart)
    options.TimeStart = 1/fs;
end

% number of EEG channels and time points
[n_chan,n_t] = size(eeg_data);


%%% SET TIME

t = (0:1:(n_t-1))/fs;
t = t+options.TimeStart;


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
% rows of Color to matcho
options.Color = rep_colors(options.Color,n_chan);

% make every other color slightly lighter or darker if ColorBeta ~= 0
if options.TintOrShade > 0
    ts_f = 1 - options.TintOrShade; % "tint/shade factor" - transform so higher values darken more
    options.Color(2:2:n_chan,:) = options.Color(2:2:n_chan,:)*ts_f;
elseif options.TintOrShade < 0
    ts_f = options.TintOrShade*-1; % negative value signifies tint; not needed for computation
    options.Color(2:2:n_chan,:) = ...
        ((1 - options.Color(2:2:n_chan,:))*ts_f)...
        + options.Color(2:2:n_chan,:); 
end


%%% OFFSET

% change offset of channels
offset_data=eeg_data;
for i=1:n_chan
    offset_data(i,:)=eeg_data(i,:)-mean(eeg_data(i,:), 'omitnan');     % first de-mean each channel so centered at zero
    offset_data(i,:)=offset_data(i,:)+options.Offset*(i-1);    % offset
end


%%% PLOT

if options.PlotNewFig
    f = figure();
    set(f,'units','centimeters','position',options.FigPos);
else
    f = [];
end

% plot
h=plot(t,offset_data','LineWidth',options.Linewidth);
set(gcf,'Color',[1 1 1]); % background colour

% change colors
for i=1:n_chan; set(h(i),'color',options.Color(i,:)); end
axis tight

% axis limits
if ~isempty(options.AxisLim)
    axis(options.AxisLim)
end

% overlay highlight if requested
if ~isempty(options.HighlightTimes)
    if size(options.HighlightTimes,2) ~= 2
        error('HighlightTimes must have two columns, with the first giving the lower bound and the second the upper bound of any highlights')
    end
    n_hi = size(options.HighlightTimes,1); % number of highlights
    
    % expand color if necessary
    options.HighlightColors = rep_colors(options.HighlightColors,n_chan);

    yl=ylim; % get y axis limits
    
    % plot
    for i=1:n_hi
        h_t = options.HighlightTimes(i,:);
        patch(h_t([1 1 2 2]),yl([1 2 2 1]),...
            options.HighlightColors(i,:),...
            'EdgeColor','none','FaceAlpha',options.HighlightAlpha)
    end
end

% plot labels and set font sizes
if ~isempty(options.ChannelNames) && ~options.ForceLabelsOff
    ax=ancestor(h(1),'axes');
    ax.YTick=0:options.Offset:(options.Offset*(n_chan-1));
    ax.YTickLabel=options.ChannelNames;
    ax.YAxis.FontSize=options.YFontSize;
    ax.XAxis.FontSize = options.XFontSize;
else
    ax=ancestor(h(1),'axes');
    ax.YTick=[];
    ax.XAxis.FontSize = options.XFontSize;
end

% label axes, title
xlabel(options.XAxisLabel,'FontSize',options.LabelFontSize);
ylabel(options.YAxisLabel,'FontSize',options.LabelFontSize);
if ~isempty(options.Title)
    title(options.Title,'FontSize',options.LabelFontSize);
end

% function for repeating colors if fewer colours than objects to colour
% will also ensure final list of colors is same size (n rows) as number of objects 
function my_colors = rep_colors(my_colors,n_final)

n_row = size(my_colors,1);
if n_row < n_final
    my_colors = repmat(my_colors,ceil(n_final/n_row),1);
end
my_colors = my_colors(1:n_final,:); 

