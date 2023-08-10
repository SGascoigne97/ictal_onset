function [cell_imprint] = imprint_to_regions(cell_imprint,json_data,data_export)
    
    for i_imprint = 1:height(cell_imprint)
        values = cell_imprint{i_imprint,1}{1};
         %%  convert imprint to regions, with threshold that I can vary easily
        channel_details = json_data(1).channel_details;
        channel_details.chan_name = strrep(channel_details.chan_name,' ','');
    
        patientChannelLabels = data_export.segment_channel_labels{1};
        patientChannelLabels = strrep(patientChannelLabels,' ','');
        no_channels = size(patientChannelLabels,1);

        for i_parcelation = 1:4
            rois = {};
            for channel = 1:no_channels
                try
                    name = patientChannelLabels(channel);
                    rois{channel} = channel_details.ROIname{strcmp(channel_details.chan_name, name)}{i_parcelation};
                catch e
                    return;
                end
            end
            
        
            rois = rois(~ismissing(rois));
            if isempty(rois)
                continue
            end
            uniqueRois = unique(rois);
            mean_imprint_region = zeros(length(uniqueRois),size(values,2));
            n_imprint_region = zeros(length(uniqueRois),1);
            for i_rois = 1:length(uniqueRois)
                mean_imprint_region(i_rois,:) = mean(values(strcmp(rois,uniqueRois(i_rois)),:),1);
                n_imprint_region(i_rois) = sum(strcmp(rois,uniqueRois(i_rois)));
            end
            warning('off', 'MATLAB:table:RowsAddedNewVars');
            cell_imprint(i_imprint, 'mean_imprint_region_' + string(i_parcelation)) = {mean_imprint_region};
            cell_imprint(i_imprint, 'n_imprint_region_' + string(i_parcelation)) = {n_imprint_region};
            cell_imprint(i_imprint, 'region_names_' + string(i_parcelation)) = {uniqueRois};
            warning('on', 'MATLAB:table:RowsAddedNewVars');
        end
    end
end

