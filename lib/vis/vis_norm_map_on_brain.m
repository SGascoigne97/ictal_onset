function f = vis_norm_map_on_brain(norm_meas,feat_names,atlas,options)
% VIS_NORM_MAP_ON_BRAIN Plot normative map (all features) on brain surface
% in one or more views.
%
% Creates a new figure with views x features subplots.
%
% Calls vis_plot_abnormalities_on_brain with options that work well for
% this application. 
%
%   f = VIS_NORM_MAP_ON_BRAIN(norm_meas,feat_names,atlas) returns the 
%   figure handle, f, to the plots of the normative map, norm_meas (ROIs x 
%   features matrix). Each column of subplots corresponds to a different 
%   feautre whose names are given in the cell array feat_names. ROI info is
%   passed via atlas, a table that must contain variables "names" (cell
%   arrays of ROI names) and "xyz" (ROI xyz coordinates). By default, the
%   first row of atlas is used for ROI info, but a different row can be
%   specified using the optional argument "AtlasIndex".
%
%   See argument validation comments for Name/Value pairs that provide 
%   limited customisation (e.g., views to plot). 
%
% See also VIS_PLOT_ABNORMALITES_ON_BRAIN.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    norm_meas {mustBeNumeric}
    feat_names cell
    atlas table
    options.View cell = {'top','anterior','left','right'};  % must be subset of the default views
    options.MarkerSize (1,1) = 100                          % size to plot ROI markers
    options.Colormap = 'parula';    % colormap color
    options.ColorbarLocation = 'southoutside';              % colorbar location; see colorbar options
    options.PlotBrain = true        % whether to plot brain surface (can be useful to just export markers in epsc format - surface needs to be exported as non-vector graphic)
    options.FontSize = 16           % fontsize (colorbar labels only)
    options.TitleFontSize = 16;     % size of title (feature names)
    options.AtlasIndex = 1;
    options.MarkerEdgeOff = true ;  % whether to turn off plotting marker edges in a different colour; LineWidthSpared will still influence total marker size
    options.LineWidth = 0.5;        % marker linewidth
    options.MarkerEdgeColor = [0 0 0]; % marker color if MarkerEdgeOff = false
end

% number of features in map
n_feat = size(norm_meas,2);

% number of views to plot brain
n_view = length(options.View);

% start figure and set size
f=figure();
set(f,'units','centimeters','position',[2 2 10+n_feat*7 10+n_view*6])

% counter for subplots
my_count=1; 

% plot for each view and feature
for i=1:n_view
    for j=1:n_feat
        
    % which rois to plot (0 = plot, NaN = missing from map)
    x = single(isnan(norm_meas(:,j)));
    x(x==1) = NaN;
        
    % only plot colorbar for last view 
    if i == n_view
        cbar = true;
    else
        cbar = false;
    end
    
    ax = subplot(n_view,n_feat,my_count);
    pos = ax.Position;
    
    % call vis_plot_abnormalities_on_brain
    vis_plot_abnormalities_on_brain(norm_meas(:,j),x,...
        atlas.names{options.AtlasIndex},atlas.xyz{options.AtlasIndex},...
        'SparedMarkerEdgeOff',true,'SparedColor',options.MarkerEdgeColor',...
        'PlotNewFig',false,'MarkerSize',options.MarkerSize,'View',options.View{i},...
        'PlotColorbar',cbar,'Colormap',options.Colormap,'ColorbarLocation',...
        options.ColorbarLocation,'FontSize',options.FontSize,'PlotBrain',options.PlotBrain,...
        'SparedMarkerEdgeOff',options.MarkerEdgeOff,'LineWidthSpared',options.LineWidth);
    
    % return subplot to full size if has colorbar
    if i == n_view
        ax.Position = pos;
    end
    
    % title
    if i == 1
        title(feat_names{j},'FontSize',options.TitleFontSize)
    end
    
    % subplot counter
    my_count = my_count+1;
    
    end
end