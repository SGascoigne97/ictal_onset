function f = vis_psd_summaries(filenames,mat_path,data_type,options)
% VIS_PSD_SUMMARIES Plots one summary PSD from each segment in either a
% single plot or a set of subplots.
%
%   f = vis_psd_summaries(filenames,mat_path,'EEG') takes a cell array of
%   file names for MATLAB workspaces containing iEEG data (filenames), all
%   located in the folder mat_path, and plots one summary PSD for each
%   file in the plot with figure handle f. 
%
%   f = vis_psd_summaries(filenames,mat_path,'PSD') takes a cell array of
%   file names for MATLAB workspaces containing PSD data (filenames), all
%   located in the folder mat_path, and plots one summary PSD for each
%   file in the plot with figure handle f. 

%   By default, the summary PSD is the mean PSD value at each frequency 
%   across all of the iEEG segment's channels.
%
%   Possible summary measures (name/value pair 'SummaryMeasure') are
%   - MeanAcross:   mean PSD value at each frequenncy across all of the segment's channels
%   - MedianAcross: median PSD value at each frequenncy across all of the segment's channels
%   - MaxAcross:    maximum PSD value at each frequenncy across all of the segment's channels
%   - MinAcross:    minimum PSD value at each frequenncy across all of the segment's channels
%
%   and
%   
%   - MeanPSD: 	 the PSD of the channel with the mean area under the PSD, relative to the segment's other channels
%   - MedianPSD: the PSD of the channel with the median area under the PSD
%   - MaxPSD:    the PSD of the channel with the maximum area under the PSD
%   - MinPSD:    the PSD of the channel with the minimum area under the PSD
%
%   Note that for the last four measures, the area under the PSD is computed
%   computed in terms of the original PSD, not in db/Hz.
%
%   Groups cab be provided as a numeric vector or cell array using the 
%   optional argument 'Groups'. The groups can either be plotted in one  
%   plot (default) or split into different subplots by group using the
%   'SplitGroups' optional argument. If plotted in the same plot
%   ('SplitGroups',false), 'Colors' and 'Labels' will be
%   used to encode and label groups rather than individual PSDs. 
%
%   See argument validation for optional name/value pair arguments. Use
%   'Interactive' to make plot interactive (does not work if 'SplitGroups'
%   is true).
%
%   See also MEAS_PWELCH_PSD, VIS_PSD, VIS_EEG
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    filenames cell
    mat_path char
    data_type char {mustBeMember(data_type,{'EEG','PSD'})}
    
    % Summary measure to compute from each set of PSDs
    options.SummaryMeasure char {mustBeMember(options.SummaryMeasure,...
        {'MeanAcross','MedianAcross','MaxAcross','MinAcross',...
        'MeanPSD','MedianPSD','MaxPSD','MinPSD'})} = 'MeanAcross'
    
    % PSD settings - only used if data_type = 'EEG'
    options.Window = 2;         % pwelch window size in seconds (see meas_pwelch_psd)
    options.Overlap = 1;        % pwelch window overlap size in seconds (see meas_pwelch_psd)
    options.Freq = 0.5:0.5:100; % frequencies at which to compute PSD
    
    % group settings
    options.Groups (:,1) = []
    options.SplitGroups = false     % whether to split groups into different subplots (vs. plotting in same plot with different colours)
    options.LinkAxes = true         % whether to link axes when groups are split; overrides 'AxisLim'
    options.NCol = 5                % maximum number of columns to use for subplots when SplitGroups is true
    
    % color
    options.Color {mustBeNumeric} = []; % colour(s) to use for lines; if groups plotted in one plot, will be used for groups. If too few, will repeat.
    options.Colormap char = 'lines';    % if options.Color is empty, name of colourmap to use
    options.ColormapN (1,1) = 7;        % if colours set using colormap, number of colours to draw 
    
    % legend/labels
    options.Labels cell = {};                           % label for each PSD; not used if Groups plotted in one plot (SplitGroups = false);
    options.NLegendCol (1,1) {mustBeNumeric} = 2;       % number of columns in legend
    options.LegendFontSize (1,1) {mustBeNumeric} = 8;   % font size for legend
    options.LegendLocation char = '';                   % if left empty, will be set in script (default is southoutside if groups that are split, and westoutside otherwise)
    options.FontSize (1,1) {mustBeNumeric} = 12;        % font size for tick labels
    options.LabelFontSize (1,1) {mustBeNumeric} = 14;   % font size for axis labels and title
    
    % other figure/plot settings
    options.Interactive {mustBeNumericOrLogical} = false;   % whether to make plot interactive (can click on and bold PSD line)
    options.PlotNewFig {mustBeNumericOrLogical} = true;     % whether plot new figure (and pass figure handle f as output); turn off to plot as subplot in existing figure
    options.FigPos = [2 2 30 20];                   % figure position in centimeters; only used for new figure; will be amended in code if SplitGroups = true
    options.Linewidth (1,1) = 1;                    % linewidth of lines used for iEEG traces
    options.AxisLim = [];                           % axis limits - should be (1,4) or axis style (e.g., 'square')

end

% number of sets of PSDs
n_pxx = length(filenames);

% initialise array for storying summary PSDs
% if PSD, load first filename to get frequency info
if isequal(data_type,'PSD')
        load([mat_path '/' filenames{1}],'freq');
        options.Freq = freq;
end
n_freq = length(options.Freq);
all_sum_pxx = zeros(n_pxx,n_freq);

%%% SUMMARY PSDs
% get summary measure for each set of PSDs
for i=1:n_pxx
    switch data_type
        % load pre-computed PSDs 
        case 'PSD'
            
            % load
            load([mat_path '/' filenames{i}],'freq','pxx');
            
            % check that frequencies are the same
            if ~isequal(freq,options.Freq)
                error('PSD frequency values do not match')
            end
        
        % load iEEG and compute PSDs
        case 'EEG'
            
            % load
            load([mat_path '/' filenames{i}],'eeg_data','eeg_fs')
            
            % compute pxx
            [pxx,~] = meas_pwelch_psd(eeg_data,...
                options.Window,options.Overlap,options.Freq,eeg_fs);

    end
    
    % summary PSD
    all_sum_pxx(i,:) = compute_summary_psd(pxx,options.SummaryMeasure);
end


%%% GROUPS

if ~isempty(options.Groups)
    % groups + number of groups
    group_names = unique(options.Groups);
    n_groups = length(group_names);
    
    % numeric encoding - will use for plotting and colours
    group_num = zeros(n_pxx,1);
    if iscell(options.Groups)
        for i=1:n_groups
            group_num(ismember(options.Groups,group_names{i})) = i;
        end
    else
        for i=1:n_groups
            group_num(options.Groups==group_names(i)) = i;
        end
    end
    
    % turn group names into strings if numeric
    if isnumeric(group_names)
        %if size(group_names,2)>1; group_names = group_names'; end % col vector
        group_names = strtrim(cellstr(num2str(group_names)));
    end
    
end



%%% COLORS
% will set and pass options.Color to vis_psd

% if no groups or groups are split into different plots, specify colour for
% each PSD; if groups in one plot, one color for each group

% number of colours needed
if ~isempty(options.Groups) && ~options.SplitGroups
    n_clr = n_groups;
else
    n_clr = n_pxx; 
end

% use Colormap if Color is empty
if isempty(options.Color)
    eval(['options.Color = ' options.Colormap '(' num2str(options.ColormapN) ');']);
    % otherwise, check that Color has three columns
elseif size(options.Color,2) ~= 3
    error('Color must either be empty or an n x 3 array')
end

% match number of PSDs or groups
options.Color = rep_colors(options.Color,n_clr);



%%% LABELS

% format title - use summary measure
plot_title = lower(options.SummaryMeasure);
plot_title = strrep(plot_title,'psd',' PSD');
plot_title = strrep(plot_title,'across',' PSD across channels');



%%% PLOT
if isempty(options.Groups) % no groups
    
    % legend location
    if isempty(options.LegendLocation)
        options.LegendLocation = 'westoutside';
    end
    
    
    f = vis_psd(all_sum_pxx,options.Freq,'PSD',...
        'Color',options.Color,...               % color
        'Interactive',options.Interactive,...   % figure settings
        'PlotNewFig',options.PlotNewFig,...
        'FigPos',options.FigPos,...
        'Linewidth',options.Linewidth,...       % linewidth
        'AxisLim',options.AxisLim,...                  % axis
        'Title',plot_title,...                   % title
        'ChannelNames',options.Labels,...       % legend info
        'NLegendCol',options.NLegendCol,...
        'LegendFontSize',options.LegendFontSize,...
        'LegendLocation',options.LegendLocation,...
        'FontSize',options.FontSize',...        % other font size info
        'LabelFontSize',options.LabelFontSize);
else

    if options.SplitGroups % groups, split into subplots
        
        % warnings for options that are ignored
        if ~options.PlotNewFig
            warning('PlotNewFig = false ignored; new figure created when Groups specified and SplitGroups is true')
        end
        if options.Interactive
            warning('Plot cannot be made interactive when Groups specified and SplitGroups is true')
        end
        
        % number of columns and rows
        n_col = min(options.NCol,n_groups);
        n_row = ceil(n_groups/n_col);
        
        % make figure
        f = figure();
        set(f,'units','centimeters','position',options.FigPos);
        
        % plot subplot for each group
        for i=1:n_groups
            
            % subplot
            ax(i) = subplot(n_row,n_col,i);
            
            % labels, handling the case of no labels provided
            if isempty(options.Labels)
                lab = {};
            else
                lab = options.Labels(group_num==i);
            end
            
            % legend location 
            if isempty(options.LegendLocation)
                options.LegendLocation = 'southoutside';
            end
            
            % plot
            vis_psd(all_sum_pxx(group_num==i,:),options.Freq,'PSD',...
                'Color',options.Color(group_num==i,:),...               % color
                'Interactive',false,...   % figure settings
                'PlotNewFig',false,...
                'Linewidth',options.Linewidth,...       % linewidth
                'AxisLim',options.AxisLim,...                  % axis
                'Title',group_names{i},...                   % title
                'ChannelNames',lab,...       % legend info
                'NLegendCol',options.NLegendCol,...
                'LegendFontSize',options.LegendFontSize,...
                'LegendLocation',options.LegendLocation,...
                'FontSize',options.FontSize',...        % other font size info
                'LabelFontSize',options.LabelFontSize);
            
        end
        sgtitle(plot_title)
        
        % link axes if requested
        if options.LinkAxes
            linkaxes(ax);
        end
        
        
    else % groups, one plot
        h_all = [];
        for i=1:n_groups
            
            plot_new = options.PlotNewFig && i==1; % if plot new fig, only do so for first group
            
            % do not plot new fig after 
            [f,h] = vis_psd(all_sum_pxx(group_num==i,:),options.Freq,'PSD',...
                'Color',options.Color(i,:),...               % color
                'Interactive',options.Interactive,...   % figure settings
                'PlotNewFig',plot_new,...  
                'FigPos',options.FigPos,...
                'Linewidth',options.Linewidth,...       % linewidth
                'AxisLim',options.AxisLim,...                  % axis
                'Title',plot_title,...                   % title
                'FontSize',options.FontSize',...        % other font size info
                'LabelFontSize',options.LabelFontSize);
            h_all(i) = h(1);
            hold on
        end
        % legend for groups
        
        % warning if labels and groups aren't split
        if ~isempty(options.Labels)
            warning('provided Labels will be overwritten by group names')
        end
        % legend location
        if isempty(options.LegendLocation)
            options.LegendLocation = 'westoutside';
        end
        
        legend(h_all,group_names,'Location',options.LegendLocation,...
            'NumColumns',options.NLegendCol,...
            'FontSize',options.LegendFontSize);
        
        hold off
    end
end

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



% compute summary PSD from set of PSDs (e.g., by averaging across channels)
function sum_pxx = compute_summary_psd(pxx,SummaryMeasure)
    
% Take mean to approximate area under each curve (probably a better way to
% do this - e.g., trap rule - but works as a quick solution)
if contains(SummaryMeasure,'PSD')
    psd_area = mean(pxx,2); % mean value of each PSD 
end

switch SummaryMeasure
    case 'MeanAcross'
        sum_pxx = mean(pxx,1);
    case 'MedianAcross'
        sum_pxx = median(pxx,1);
    case 'MaxAcross'
        sum_pxx = max(pxx,[],1);
    case 'MinAcross'
        sum_pxx = min(pxx,[],1);
    case 'MeanPSD'
        % find psd with area closest to mean area 
        mean_area = mean(psd_area);
        area_diff = psd_area - mean_area;
        [~,idx] = min(abs(area_diff));
        sum_pxx = pxx(idx,:);
    case 'MedianPSD'
         % find psd with area closest to median area 
        med_area = median(psd_area);
        area_diff = psd_area - med_area;
        [~,idx] = min(abs(area_diff));
        sum_pxx = pxx(idx,:);
    case 'MaxPSD'
        [~,idx] = max(psd_area);
        sum_pxx = pxx(idx,:);
    case 'MinPSD'
        [~,idx] = min(psd_area);
        sum_pxx = pxx(idx,:);
end
end

