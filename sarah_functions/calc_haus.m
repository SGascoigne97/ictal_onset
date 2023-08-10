function [comp_table] = calc_haus(pat_onset, atlas_scale, opts)
% Compute average EEG across regions of interest for each patient

% Sarah J Gascoigne 10/02/2023

% input:
%   - pat_onset: onset_output table for one specific patient 
%   - atlas_scale: table containing ROI atlas at chosen scale
%   - optional inputs
%       - det_method: state the method of onset detection you'd like to use
%           ("imprint", "EI", "PLHG")
%       - comparison: state whether the comparisons are based on
%           comparisons between onset and resected regions ("resction") or
%           between pairs of seizure onsets ("pairwise")

% output
%   - comp_table: table containing raw and normalised Hausdorff's distance

    arguments
        pat_onset
        atlas_scale 
        opts.det_method (1,1) string {mustBeMember(opts.det_method, ["CLO", "imprint", "EI", "PLHG"])} = "imprint" 
        opts.comparison (1,1) string {mustBeMember(opts.comparison, ["resection", "pairwise"])} = "resection" 
        opts.chan_or_roi (1,1) string {mustBeMember(opts.chan_or_roi, ["chan", "roi_120", "roi_250"])} = "roi_120" 
    end
    
    %fill in optional arguments
    det_method = opts.det_method;
    comparison = opts.comparison;
    chan_or_roi = opts.chan_or_roi;

    % Pull in xyz coordinates from Atlas
    xyz = atlas_scale.xyz{1,1};

    % Pull out the unique ROIs recorded for this patient
    unq_roi = pat_onset.(sprintf("roi_names_%s", extractAfter(chan_or_roi, "roi_"))){1,1}; % Need to find out if I can get distances between channels (across grey matter) 

    if det_method == "CLO"
          % Extract the automatically detected onsets matrix
        onset_binary = cell2mat(pat_onset.sprintf("clo_%s", chan_or_roi));
        seizure_ids = "CLO";
        sz_count = 1;
    else
        % Extract the automatically detected onsets matrix
        onset_binary = cell2mat(pat_onset.(sprintf('%s_%s',det_method, chan_or_roi)));
        % Extract seizure IDs
        seizure_ids = string(pat_onset.Segment_ids{1,1});
        
        % Remove seizures with no onset detected
        % This should have been done earlier, double check code!
        %onset_binary = onset_binary(:,sum(onset_binary) > 0);
        %seizure_ids = seizure_ids(sum(onset_binary) > 0);
    
        sz_count = size(onset_binary,2);
    end

    if comparison == "resection"
        
        comp_table = array2table(zeros(sz_count,1), 'VariableNames', ...
            {'Hausdorff'});
% 
%         % Obtain xyz coordinated for resected regions
%         resec_binary = ismember(strrep(string(atlas_scale.name{1,1}),' ',''),...
%             strrep(unq_roi(find(resected_binary)), ' ',''));
%         nan_resec = isnan(resected_binary); 
%         resec_xyz= xyz(resec_binary,:);

        % Iterate through seizures and compute Hausdorff's distance
        for sz = 1:sz_count
            % Add in patient ID
            comp_table.Patient_id(sz) = pat_onset.Patient_id;
            % Add in seizure ID
            comp_table.Seizure_id(sz) = seizure_ids(sz);
            
            % Obtain onset for this seizure
            sz_onset= onset_binary(:,sz);

            % Obtain xyz coordinated for resected regions
            % Extract resected region for this patient
            resected_binary = cell2mat(pat_onset.(sprintf('resected_%s', chan_or_roi)));
            nan_resec = isnan(resected_binary); 

            % If any regions are ignored, remove them from analysis now 
            loc_nan = isnan(sz_onset)+nan_resec >0;
            sz_onset_rm_nan = sz_onset(~loc_nan);
            unq_roi_rm_nan = unq_roi(~loc_nan);
            resected_binary_rm_nan = resected_binary(~loc_nan);
            % locate regions in atlas
            resec_binary = ismember(strrep(string(atlas_scale.name{1,1}),' ',''),...
            strrep(unq_roi_rm_nan(find(resected_binary_rm_nan)), ' ',''));
            resec_xyz= xyz(resec_binary,:);

            % Ignore seizures with no onset detected
            if sum(sz_onset_rm_nan, 'omitnan')==0 || isempty(resec_xyz)
                comp_table.Hausdorff(sz) = NaN;
                comp_table.Hausdorff_norm(sz) = NaN;
                continue
            end

            % Compute Hausdorff distance between onset regions and resection
            onset_roi_binary = ismember(strrep(string(atlas_scale.name{1,1}),' ',''),...
                strrep(unq_roi_rm_nan(find(sz_onset_rm_nan)),' ',''));
           
            onset_xyz = xyz(onset_roi_binary,:);
           
            [hd,~] = HausdorffDist(resec_xyz, onset_xyz);
      
            comp_table.Hausdorff(sz) = hd;
            % Scale Hausdorff's distance to between zero and one by dividing by the
            % maximum distance between two ROIs in the atlas
            comp_table.Hausdorff_norm(sz) = hd/max(max(atlas_scale.dists{1,1}));
            
        end


        % Reorder table to help with readability
        comp_table = comp_table(:,{'Patient_id','Seizure_id', 'Hausdorff', 'Hausdorff_norm'});
        
    elseif comparison == "pairwise"
        if det_method == "CLO"
            print("Pairwise comparison is not possible for clinically labelled onsets")
            comp_table = [];
        else
            comp_table = array2table(nan(sum(1:(sz_count-1)),1), 'VariableNames', ...
                {'Hausdorff'});
            count = 0;
            for sz1 = 1:sz_count
                for sz2 = 1:sz_count
                    if sz1<sz2
                        count = count+1;
                         comp_table.Seizure1_id(count) = seizure_ids(sz1);
                         comp_table.Seizure2_id(count) = seizure_ids(sz2);
                         % Compute Hausdorff distance between onset regions for
                         % seizures one and two

                         % If either seizure has no onset detected, skip
                         % this comparison
                         sz_onset1 = onset_binary(:,sz1);
                         sz_onset2 = onset_binary(:,sz2);
                            % Remove any ignored regions
                         na_rm = isnan(sz_onset1) + isnan(sz_onset2);

                         sz_onset1 = sz_onset1(~na_rm);
                         sz_onset2 = sz_onset2(~na_rm);
                         unq_roi_rm_nan = unq_roi(~na_rm);

                        if sum(sz_onset1, 'omitnan')==0 || sum(sz_onset2, 'omitnan')==0 
                            comp_table.Hausdorff(count) = NaN;
                            comp_table.Hausdorff_norm(count) = NaN;  
                            % Add in patient ID
                            comp_table.Patient_id(count) = pat_onset.Patient_id;
                            continue
                        end

                         onset1_roi_binary = contains(string(atlas_scale.name{1,1}), unq_roi_rm_nan(find(sz_onset1)));
                         onset1_xyz = xyz(onset1_roi_binary,:);
                         onset2_roi_binary = contains(string(atlas_scale.name{1,1}), unq_roi_rm_nan(find(sz_onset2)));
                         onset2_xyz = xyz(onset2_roi_binary,:);
                         [hd,~] = HausdorffDist(onset1_xyz, onset2_xyz);
                         comp_table.Hausdorff(count) = hd;
                         % Scale Hausdorff's distance to between zero and one by dividing by the
                         % maximum distance between two ROIs in the atlas
                         comp_table.Hausdorff_norm(count) = hd/max(max(atlas_scale.dists{1,1}));
                         % Add in patient ID
                        comp_table.Patient_id(count) = pat_onset.Patient_id;
    
                    end
                end
            end
    
        % Reorder table to help with readability
        comp_table = comp_table(:,{'Patient_id', 'Seizure1_id', 'Seizure2_id', 'Hausdorff', 'Hausdorff_norm'});
        end
        
    end
end
