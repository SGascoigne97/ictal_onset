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
    end
    
    %fill in optional arguments
    det_method = opts.det_method;
    comparison = opts.comparison;

    % Pull in xyz coordinates from Atlas
    xyz = atlas_scale.xyz{1,1};

    % Pull out the unique ROIs recorded for this patient
    unq_roi = pat_onset.ROI_ids{1,1};

    if det_method == "CLO"
          % Extract the automatically detected onsets matrix
        onset_binary = cell2mat(pat_onset.Labelled_onset);
        seizure_ids = "CLO";
        sz_count = 1;
    else
        % Extract the automatically detected onsets matrix
        onset_binary = cell2mat(pat_onset.(sprintf('%s_roi',det_method)));
        % Extract seizure IDs
        seizure_ids = string(pat_onset.Segment_ids{1,1});
        
        % Remove seizures with no onset detected
        % This should have been done earlier, double check code!
        onset_binary = onset_binary(:,sum(onset_binary) > 0);
        seizure_ids = seizure_ids(sum(onset_binary) > 0);
    
        sz_count = size(onset_binary,2);
    end

    if comparison == "resection"
        % Extract resected region for this patient
        resected_binary = cell2mat(pat_onset.Resected);
        comp_table = array2table(zeros(sz_count,1), 'VariableNames', ...
            {'Hausdorff'});

        % Obtain xyz coordinated for resected regions
        resec_binary = ismember(strrep(string(atlas_scale.name{1,1}),' ',''),...
            strrep(unq_roi(find(resected_binary)), ' ',''));
        resec_xyz = xyz(resec_binary,:);

        % Iterate through seizures and compute Hausdorff's distance
        for sz = 1:sz_count
            % Add in patient ID
            comp_table.Patient_id(sz) = pat_onset.Patient_id;
            % Add in seizure ID
            comp_table.Seizure_id(sz) = seizure_ids(sz);
            % Compute Hausdorff distance between onset regions and resection
            onset_roi_binary = ismember(strrep(string(atlas_scale.name{1,1}),' ',''),...
                strrep(unq_roi(find(onset_binary(:,sz))),' ',''));
            if sum(onset_roi_binary)==0
                comp_table.Hausdorff(sz) = NaN;
                comp_table.Hausdorff_norm(sz) = NaN;
                continue
            end
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
                         onset1_roi_binary = contains(string(atlas_scale.name{1,1}), unq_roi(find(onset_binary(:,sz1))));
                         onset1_xyz = xyz(onset1_roi_binary,:);
                         onset2_roi_binary = contains(string(atlas_scale.name{1,1}), unq_roi(find(onset_binary(:,sz2))));
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
