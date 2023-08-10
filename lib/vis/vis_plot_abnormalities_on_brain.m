function f = vis_plot_abnormalities_on_brain(z,roi_is_resect,roi_names,roi_xyz,options)
% VIS_PLOT_ABNORMALITIES_ON_BRAIN Plot ROI abnormalities on brain surface, 
% viewed from top, left, or right, or anterior. For left and right views, 
% will only plot left/right ROIs.
%
% By default, resected ROIs will be outlined in black and spared ROIs in
% white.
% 
%   f = VIS_PLOT_ABNORMALITIES_ON_BRAIN(z,roi_is_resect,roi_names,roi_xyz)
% 
%   See argument validation comments for variable info and name/value pair
%   arguments.
% 
% Requires fieldtrip to load brain surface.
% 
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% October 2022; modified January 2023

arguments
    z {mustBeNumeric}   % vector of ROIs abnonormalities to plot, length n ROIs
    roi_is_resect       % boolean of which ROIs were resected, length n
    roi_names           % cell array of names of ROIs in z, length n - used to determine L (contains "l." or "Left-") and R (contains "r." or "Right-") hemisphere ROIs. 
    roi_xyz             % ROI xyz coordinates, size n x 3
    options.MarkerSize (1,1) = 300  % size to plot ROI markers
    options.CLim = [min(z) max(z)]  % colorbar limits
    options.PlotBrain = true        % whether to plot brain surface (can be useful to just export markers in epsc format - surface needs to be exported as non-vector graphic)
    options.View {mustBeMember(options.View,{'top','left','right','anterior'})} = 'top' % view of brain; left and right will only plot ROIs in that hemisphere
    options.LineWidthSpared {mustBeGreaterThan(options.LineWidthSpared,0)} = 1.5       % linewidth of outside of spared ROI markers
    options.LineWidthResected {mustBeGreaterThan(options.LineWidthResected,0)} = 1.5;  % linewidth of outside of resected ROI markers
    options.FontSize = 18               % fontsize (colorbar labels only)
    options.SparedColor = [1 1 1];     % colour of outine of spared ROIs
    options.ResectedColor = [0 0 0];   % colour of outline of resected ROIs
    options.PlotNewFig = true;          % whether plot new figure (and pass figure handle f as output) - can turn off if, e.g., want to plot in subplots of existing figure
    options.PlotColorbar = true;        % whether to plot colorbar
    options.surface_path char = 'lib/vis/brain_surf';   % folder containing lh.pial and rh.pial
    options.all_xyz = roi_xyz           % coordinates of ALL ROIs in parcellation, used to set axis limits, if different from roi_xyz
    options.SparedMarkerEdgeOff = false;% whether to turn off plotting marker edges in a different colour on spared ROIs; LineWidthSpared will still influence total marker size
    options.ColorbarLocation = 'eastoutside'; % colorbar location; see colorbar options
    options.Colormap = 'parula'; % colormap color
end

% load surface 
surface_path = options.surface_path;
[lhpialvert,lhpialface]=read_surf([surface_path '/lh.pial']);
lhpialface=lhpialface+1;
[rhpialvert,rhpialface]=read_surf([surface_path '/rh.pial']);
rhpialface=rhpialface+1;
lhpialvert(:,2)=lhpialvert(:,2)-18;
lhpialvert(:,3)=lhpialvert(:,3)+20;
rhpialvert(:,2)=rhpialvert(:,2)-18;
rhpialvert(:,3)=rhpialvert(:,3)+20;

% mark left and right hemisphere ROIs
left_roi = contains(roi_names,'Left-') | contains(roi_names,'l.');
right_roi = contains(roi_names,'Right-') | contains(roi_names,'r.');

% only keep left/right ROIS for left/right views
n_roi = length(z);
switch options.View
    case 'left'
        keep_roi = left_roi;
    case 'right'
        keep_roi = right_roi;
    otherwise
        keep_roi = true(1,n_roi);
end
z = z(keep_roi);
roi_is_resect = roi_is_resect(keep_roi);
roi_xyz = roi_xyz(keep_roi,:);

% how much to shift x-axis depending on view
switch options.View
    case 'top'
        roi_shift_x = 0;
        roi_shift_y = 0;
        roi_shift_z = 140;
    case 'left'
        roi_shift_x = -140;
        roi_shift_y = 0;
        roi_shift_z = 0;
    case 'right'
        roi_shift_x = 140;
        roi_shift_y = 0;
        roi_shift_z = 0;
    case 'anterior'
        roi_shift_x = 0;
        roi_shift_y = 140;
        roi_shift_z = 0;       
end

% marker size
mk_sz = options.MarkerSize;

% resected/spared booleans - will exclude any ROI with roi_is_resect = NaN
resected = roi_is_resect==1;
spared = roi_is_resect==0;

% plot
if options.PlotNewFig
    f=figure();
else
    f=[];
end
set(gcf,'Color','w')
hold on
% plot spared
if options.SparedMarkerEdgeOff % no marker edges
    scatter3(roi_xyz(spared,1)+roi_shift_x,...
        roi_xyz(spared,2)+roi_shift_y,...
        roi_xyz(spared,3)+roi_shift_z,...
        mk_sz,z(spared),'filled','LineWidth',options.LineWidthSpared);
else
scatter3(roi_xyz(spared,1)+roi_shift_x,...
        roi_xyz(spared,2)+roi_shift_y,...
        roi_xyz(spared,3)+roi_shift_z,...
        mk_sz,z(spared),'filled','MarkerEdgeColor',options.SparedColor,...
        'LineWidth',options.LineWidthSpared);
end

% plot resected
scatter3(roi_xyz(resected,1)+roi_shift_x,...
    roi_xyz(resected,2)+roi_shift_y,...
    roi_xyz(resected,3)+roi_shift_z,...
    mk_sz,z(resected),'filled','MarkerEdgeColor',options.ResectedColor,...
    'LineWidth',options.LineWidthResected);

% colour settings
colormap(options.Colormap)
if options.PlotColorbar
    colorbar(options.ColorbarLocation)
end
caxis(options.CLim)

% plot view and format
axis equal; grid off; axis off
all_xyz = options.all_xyz;
switch options.View
    case 'top'
        view([0 90]);
        xlim([min(all_xyz(:))-20 max(all_xyz(:))+20]);
        ylim([min(all_xyz(:))-20 max(all_xyz(:))+20]);
        zlim([min(all_xyz(:))-20 max(all_xyz(:))+140]);
    case 'left'
         xlim([min(all_xyz(:,1))-140 max(all_xyz(:,1))+20]);
         ylim([min(all_xyz(:,2))-20 max(all_xyz(:,2))+20]);
         zlim([min(all_xyz(:,3))-20 max(all_xyz(:,3))+20]);
        view([270 0]);
    case 'right'
        xlim([min(all_xyz(:,1))-20 max(all_xyz(:,1))+140]);
        ylim([min(all_xyz(:,2))-20 max(all_xyz(:,2))+20]);
        zlim([min(all_xyz(:,3))-20 max(all_xyz(:,3))+20]);
        view([90 0]);
    case 'anterior'
        view([180 0]);
        xlim([min(all_xyz(:))-5 max(all_xyz(:))+15]);
        ylim([min(all_xyz(:))-20 max(all_xyz(:))+140]);
        zlim([min(all_xyz(:))-20 max(all_xyz(:))+20]); 
        
end
set(gca,'FontSize',options.FontSize)

% plot brain surface
if options.PlotBrain
    switch options.View
        case 'left'
            trisurf(lhpialface,lhpialvert(:,1),lhpialvert(:,2),lhpialvert(:,3),'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)
        case 'right'
            trisurf(rhpialface,rhpialvert(:,1),rhpialvert(:,2),rhpialvert(:,3),'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)
        case {'top', 'anterior'}
            trisurf(lhpialface,lhpialvert(:,1),lhpialvert(:,2),lhpialvert(:,3),'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)
            trisurf(rhpialface,rhpialvert(:,1),rhpialvert(:,2),rhpialvert(:,3),'facecolor',[0.8 0.8 0.8],'facealpha',0.2,'edgealpha',0)
            
    end
end

% plot L/R for anterior and top view
ft = 20;
switch options.View
    case 'anterior'
        text(-40,-20,90,'L','FontSize',ft)
        text(60,-20,90,'R','FontSize',ft)
    case 'top'
        text(-80,50,90,'L','FontSize',ft)
        text(60,50,90,'R','FontSize',ft)

end

