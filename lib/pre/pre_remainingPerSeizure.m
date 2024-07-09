function [sld, thresholdSetCheck] = pre_remainingPerSeizure(pat_json, remaining_per_sz)
    % Creates a figure showing each seizure by how many noisy channels they have
    % Allows the user to drag a slider setting a threshold value that can
    % then be used to remove seizures that have a lot of noisy channels
    fprintf('Patient %s:\n',...
        pat_json(1).readable_patID)
    
    % Visualise seizure index against number of noisy channels
    % figure size made using trial and error, might need changing when
    % tested on other screens
    figure('units','pixels','Position',[2151 141 941 700]);
    colormap(jet(4))
    scatter(remaining_per_sz, 1:length(pat_json), 40,'filled');
    xlim([(min(remaining_per_sz) - 5), (max(remaining_per_sz)+1)])
    ylabel('Seizure index')
    xlabel('Count of remaining channels')
    title(sprintf("Patient %s number of remaining channels per seizure", pat_json(1).readable_patID))
    % Create slider at bottom of graph, offsetting from the axes position
    sld = uicontrol('Style', 'slider',...
        'Min',(min(remaining_per_sz) - 5),'Max',max(remaining_per_sz)+1,'Value',(min(remaining_per_sz) - 5),...
        'Units', 'Normalized',...
        'Position', gca().Position + [-0.02 -0.87 0.04 0],...
        'Callback', @thresholdChanged,...
        'SliderStep', [1/(max(remaining_per_sz)+1) 0.01]);
    % Create checkbox next to slider
    thresholdSetCheck = uicontrol('Style','checkbox','String','Threshold set', 'Position',[10 22 95 20]);
    
        function thresholdChanged(hObject,callbackdata)
        % function to force value to integer and redraw plot highlighting
        % seizures to be removed
        newval = hObject.Value;                         %get value from the slider
        newval = round(newval);                         %round off this value
        set(hObject, 'Value', newval);                  %set slider position to rounded off value
        % redraw scatter plot highlighting seizures to be removed in red
        if newval <= max(remaining_per_sz)
            scatter(remaining_per_sz, 1:length(pat_json), 40, (remaining_per_sz <= newval) + 1, 'filled');
        else
            scatter(remaining_per_sz, 1:length(pat_json), 40,'filled');
        end
        xlim([(min(remaining_per_sz) - 5), (max(remaining_per_sz)+1)])
        ylabel('Seizure index')
        xlabel('Count of remaining channels')
        title(sprintf("Patient %s number of remaining channels per seizure", pat_json(1).readable_patID))
        disp(['Threshold set to ' num2str(newval)]);     %display the value pointed by slider
    end
end