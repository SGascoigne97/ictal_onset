function [data_tbl] = av_across_roi(data_tbl, opts)
% Compute average EEG across regions of interest for each patient

% input:
%   - data_tbl: full data table
%   - optional inputs
%       - sz_type: seizure type of interest
%       - min_sz_duration: minimum duration for each seizure 
%       - min_sz_count: minimum number of seizures to ahve been recorded per patient

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl
        opts.av_method (1,1) string {mustBeMember(opts.av_method, ["median", "mean"])} = "median" % declare how EEG should be averaged (median or mean)
        opts.show_fig (1,1) logical = 1 % 1 if figure should be created
        opts.vis_sz (1,1) double {mustBeNumeric} = 10; % Choose which seizure to present (as a index in the data table)
    end
    
    %fill in optional arguments
    av_method = opts.av_method;
    show_fig = opts.show_fig;
    vis_sz = opts.vis_sz;

    % Plot channel-wise data if optional argument asks for plot
    if show_fig == 1
        sz_data = data_tbl(vis_sz,:);
        figure()
        offset = 500;
        subplot(2,1,1)
        for r = 1:size(sz_data.segment_data{1,1},1)
            plot((1:size(sz_data.segment_data{1,1},2))/256,sz_data.segment_data{1,1}(r,:)+offset)
            xlim([0, (240+sz_data.duration)])
            hold on
            offset = offset + 500;
        end
        title('Channels')
    end

    %% Average EEG time series across ROI 
    % Iterate through seizures 
    for sz = 1:size(data_tbl, 1)
        % Remove missing ROIs
        incl_roi = cellfun(@ischar,data_tbl(sz,:).channel_roi{1,1});
        roi_data = data_tbl(sz,:).channel_roi{1,1};
        rois = string(roi_data(incl_roi));
        unq_roi = unique(rois);
        roi_size = zeros(length(unq_roi),1);
        roi_data = zeros(length(unq_roi), size(data_tbl(sz,:).segment_data{1,1},2));
        for grp = 1:length(unq_roi)
            if av_method == "median"
                roi_data(grp,:) = median(data_tbl(sz,:).segment_data{1,1}(rois == unq_roi(grp),:),1);
            elseif av_method == "mean"
                roi_data(grp,:) = mean(data_tbl(sz,:).segment_data{1,1}(rois == unq_roi(grp),:),1);
            end
            roi_size(grp) = size(data_tbl(sz,:).segment_data{1,1}(rois == unq_roi(grp),:),1);
        end
        data_tbl.segment_data{sz,1} =  roi_data;
    end
    

    % If the optional argument to create a figure is 1
    if show_fig == 1
        offset = 500;
        sz_data = data_tbl(vis_sz,:);
        subplot(2,1,2)
        for r = 1:size(sz_data.segment_data{1,1},1)
            plot((1:size(sz_data.segment_data{1,1},2))/256,sz_data.segment_data{1,1}(r,:)+offset)
            xlim([0, (240+sz_data.duration)])
            hold on
            offset = offset + 500;
        end
        title('ROIs')
        sgtitle(sprintf('Patient %s (Sz %s): Comparing channel activity to medain ROI activity', sz_data.patient_id, string(sz_data.segment_id)))
    end
end
