function  [tbl_imprint_out, cell_imprint, cell_t, cell_madscores,cell_pre_features_mad, cell_pre_mahal_mat]  =...
mahal_imprint(data_tbl,sz_mat_tab,t,opts)

% calculates "imprint" and other (initial) recruitment metrics of seizures.
% 
% inputs:
%   - data_tbl: table with columns for ids,pre,duration,fs, and data for better
%       plots
%   - val_cell: cell array with data matrix chns-by-time for each segment,
%       cell columns are different features, rows are different segments
%   - t: cell array with time arrays for each segment where time
%       arrays contain time points in seconds corresponding to the columns of
%       data matrix of that segment
%   - opts: for more info on parameters look below
%outputs are in order of data_tbl, as we assume that input is in order of
%data_tbl

    arguments
        data_tbl
        sz_mat_tab
        t
        opts.rec_thresh (1,1) double {mustBeNumeric} = 9
        opts.ict_buffer (1,1) double {mustBeInteger, mustBePositive} = 10 %in seconds
        opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
        opts.prop_rec (1,1) double {mustBeNumeric, mustBePositive} = 0.7 %
    end
    
    %fill in optional arguments
%     rec_type = opts.rec_type;
    rec_thresh=opts.rec_thresh;    
    mad_thresh=opts.mad_thresh;
    ict_buffer=opts.ict_buffer;
    prop_rec = opts.prop_rec;
    
    mc=-1/(sqrt(2)*erfcinv(3/2)); % fixed factor for MAD score calculation
    %% fill in some info from data_tbl

%     nch = length(data_tbl.segment_channel_labels{1});%!!!!!!!
%     fs=data_tbl.segment_fs(1);
%     nsz = size(data_tbl,1);%number of seizures
    
    if size(data_tbl,1) ~= size(sz_mat_tab,1) || size(data_tbl,1) ~= numel(t) || ...
        size(sz_mat_tab,1) ~= numel(t)
        error('All inputs must have same amount of rows.')
    end
    
    nsegs = size(data_tbl,1);
    nmaxchr = zeros(nsegs,1);
    tmaxchr = zeros(nsegs,1);
    cell_infos = cell(nsegs,1);
    %nmaxchr = zeros(nsegs,1); % Commented out for now as this is not used in the function
    %tmaxchr = zeros(nsegs,1); % Commented out for now as this is not used in the function
    cell_madscores = table(cell(nsegs,1), cell(nsegs,1),...
        cell(nsegs,1), cell(nsegs,1), 'VariableNames', ["mahal_MAD",...
        "Forward sum", "Backward sum", "Middle sum"]);
    cell_pre_mahal_mat = cell(nsegs,1);
    cell_t = cell(nsegs,1);
    cell_pre_features_mad = cell(nsegs,1);

    patient = sprintf("UCLH%s", string(data_tbl(1,:).patient_id));

    %% Create empty output tables
   tbl_imprint_out = table(repelem(patient,nsegs,1), data_tbl.segment_id,...
        nan(nsegs,1), cell(nsegs,1), nan(nsegs,1),...
        'VariableNames',["patient_id", "segment_id", "Onset time", "Onset",...
        "num_chan_imprint"]);

   cell_imprint = table(data_tbl.segment_id, cell(nsegs,1), 'VariableNames',...
       ["segment_id", "cell_imprint"]);

    %% Iterate through seizures 
    for s = 1:nsegs
        fprintf("\n Seizure %d", s)
        %read out seizure segment & compress all features into one 3D array
       
        tw=t{s};%timing for each seizure, in actual seconds
        pre_ids = find(tw<=-(1*ict_buffer+1));
        pre_vals = sz_mat_tab.feat_mat{s}(:,pre_ids,:);
        ictal_ids = find(tw>-1*ict_buffer& tw<=data_tbl.duration(s));
        ict_vals = sz_mat_tab.feat_mat{s}(:,ictal_ids,:);
        
%         %% Remove preictal outliers (removed here in favour of removing
%         based on MAD, see lines 132-137)
%         for chan = 1:size(pre_vals,1)
%             for mark = 1:8
%                 % Remove preictal outliers (e.g. spikes)
%                 pre_vals(chan,:,mark) = filloutliers(pre_vals(chan,:,mark), NaN);
%             end
%         end
%         
        %% For each channel create a distribution of Mahalanobis distances (remove
        % time point and compare distance against all other time points)
        % Create a distribution of Mahalanobis distances (one value per
        % channel, marker, and time point)
        pre_mahal_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        for chan = 1:size(pre_vals,1)
            pre_mahal_mat(chan, :) = mahal(squeeze(pre_vals(chan,:,:)), squeeze(pre_vals(chan,:,:)));
            % For efficiency, we are not using a leave-one-out approach
            % here, the code below will perform this but makes the function
            % extremely slow

            % pre_time_points = 1:size(pre_vals(chan,:,1),2);
            % for pre_time_point = pre_time_points
            %     if all(~isnan(squeeze(pre_vals(chan,pre_time_point,:))))
            %         ref_dist = squeeze(pre_vals(chan,pre_time_points(pre_time_points~=pre_time_point),:));
            %         ref_dist = ref_dist(~isnan(sum(ref_dist,2)),:);
            %         pre_mahal_mat(chan, pre_time_point) = mahal(squeeze(pre_vals(chan,pre_time_point,:))', ref_dist);
            %     end
            % end
        end
        
          %% MAD score ictal Mahalanobis distances against preictal Mahalanobis distances
        ict_mahal_mat = nan(size(ict_vals,1), size(ict_vals(1,:,1),2));
        pre_mahal_mad_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        mahal_mad_mat = nan(size(ict_vals,1), size(ict_vals(1,:,1),2));
        
        for chan = 1:size(pre_vals,1) 
            pre_dist = pre_mahal_mat(chan,:);
            % Compute preictal median Mpre_mahal_matahalanobis distance
            pre_med = median(pre_dist, 2, 'omitnan');
            pre_smad=mc*mad_rewrite(pre_dist,1,2);
            % Compute preictal MAD values
            pre_mahal_mad_chan = (pre_dist-pre_med)./pre_smad;
            % Remove any preictal windows with MAD >= MAD threshold
            % (essentially removing preictal noise)
            ref_vals = squeeze(pre_vals(chan,:,:));
            ref_vals(pre_mahal_mad_chan >= mad_thresh,:) = [];
            pre_mahal_mad_chan(pre_mahal_mad_chan >= mad_thresh) = NaN;

            pre_mahal_mad_mat(chan,:) = pre_mahal_mad_chan;%score preictal to median & scaled mad
            % Remove preictal outliers from pre_dist before computing MAD
            pre_dist = mahal(ref_vals, ref_vals);
            % Recompute med and smad
            pre_med = median(pre_dist, 1, 'omitnan');
            pre_smad=mc*mad_rewrite(pre_dist',1,2);
            % store updated pre_mahal_mat for export 
            pre_mahal_mat(chan,~isnan(pre_mahal_mad_chan)) = pre_dist;
            pre_mahal_mat(chan,isnan(pre_mahal_mad_chan)) = NaN;
            % Compute ictal MAD values (ictal values MAD scored against
            % preictal distribution)
            ict_mahal_mat(chan, :) = mahal(squeeze(ict_vals(chan,:,:)), ref_vals);
            mahal_mad_mat(chan,:) = (ict_mahal_mat(chan,:)-pre_med)./pre_smad; %score ictal to median & scaled mad
        end

        
        %%
%         figure()
%         imagesc(mahal_mad_mat>=5)
%         %
%         markers = ["line length", "energy", "delta", "theta", "alpha", "beta", "low gamma", "high gamma"];
%         chan = 23;
%         figure()
%         tiledlayout(8,8)
%         for mark1 = 1:8
%             for mark2 = 1:8
%                 if mark1 < mark2
%                    nexttile
%                    scatter(pre_vals(chan,:,mark1), pre_vals(chan,:,mark2),10,'.')
%                    hold on
%                    scatter(ict_vals(chan,1:10,mark1), ict_vals(chan,1:10,mark2),10,'.')
%                    xlabel(markers(mark1))
%                    ylabel(markers(mark2))
%                    hold off
%                 elseif mark1 == mark2
%                     nexttile
%                     axis off
%                 else
%                     nexttile
%                     axis off
%                 end
%         
%             end
%         end
        
        %%
    %     time_points = (0:9)/8;
    %     figure("Position", [0,0,1050,1500])
    %     tiledlayout(5,1)
    %     for chan = 61:65
    %         nexttile
    %         histogram(pre_mahal_mat(chan,:))
    %         title(string(data_tbl(1,:).segment_channel_labels{:}(chan)))
    %         xlabel("Preictal Mahalanobis distance")
    %         ylabel("Frequency")
    %         xlim([0, 45])
    %         for ict_time_point = 1:10
    %             xline(ict_mahal_mat(chan,ict_time_point), 'r',...
    %                 sprintf("t(%d), MAD = %.3f",...
    %                 ict_time_point, mahal_mad_mat(chan, ict_time_point)))
    %         end
    %     end
    %   
        
        %% Add inclusion criteria
        wl=(tw(2)-tw(1));
        recruitment_threshold = rec_thresh/wl;
       
        ms_a=movsum(mahal_mad_mat>=mad_thresh,[0 recruitment_threshold-1],2);% forward looking sum
        ms_b=movsum(mahal_mad_mat>=mad_thresh,[recruitment_threshold-1 0],2);% backward looking sum
        ms_c=movsum(mahal_mad_mat>=mad_thresh,[(recruitment_threshold/2)-1 (recruitment_threshold/2)],2); % sum across centre
        imprint = ms_a >= recruitment_threshold*prop_rec | ms_b >= recruitment_threshold*prop_rec | ms_c >= recruitment_threshold*prop_rec ;%combining the backward and forward looking sum has the advantage of not cutting off early activity, or neglecting late activity
        
        nchr=sum(sum(imprint,2)>=1);%get number of channels ever in imprint
        onset_time = find(sum(imprint),1, 'first');
        if isempty(onset_time)
            onset = [];
            onset_time = NaN;
        else
            onset = sum(imprint(:,onset_time+0:7),2)>=1;
        end

        % Pull out ictal segment
        % ict_features= features(:,ictal_ids,:);
        % ict_features_mad = abs(mad_score_features(:,ictal_ids,:));

        % write output
        cell_imprint(s,:) = table(data_tbl.segment_id(s), {imprint});
        cell_infos{s}={pre_mahal_mat,mahal_mad_mat,ms_a,ms_b,ms_c};
        % cell_pre_features_mad{s} = feat_pre_smad;
        cell_madscores(s,:) = table({mahal_mad_mat},...
            {ms_a}, {ms_b}, {ms_c});
        cell_pre_mahal_mat{s} = pre_mahal_mat;
        cell_pre_features_mad{s} = pre_mahal_mad_mat;
        cell_t{s}=tw(ictal_ids);

        tbl_imprint_out(s,3:5) = table((onset_time-1)/8, {onset}, nchr);
        
    
    end
%tbl_imprint_out=table(nchr, nmaxchr, tmaxchr,'VariableNames',{'num_chan_imprint','num_chan_max_concurrent','time_till_max_concurrent'});
    
end
