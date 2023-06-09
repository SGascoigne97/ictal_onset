 

% Create table with comparison measures 

% Parameters
%   - Atlas (Lausanne 120 and 250)
%   - Summarise seizures? (Per seizure or summarised across seizures)
%   - Detection method (Imprint, EI, PLHG) 
%   - Seizure types (focal and/or subclinical) - Seozure type is not
%   included yet, this will be added later
%   - Proportion of seizures to included in across-seizures summary (0.01
%   if just one seizure is required) 
%       Note that if any other value is used, we will need to exclude
%       patients with too few seizures (e.g., 0.25 - only patients with >=4
%       seizures will be kept)

% Parameters to be scanned to determine if results are robust
final_output_both_atlas = final_output;
atlas = [120, 250];
per_sz_or_all_sz = ["all_sz", "per_sz"];
det_method = ["imprint", "EI", "PLHG"];

% Paramaters to select at start (will not be scanned)
sz_types = "all"; %["focal", "subclin", "all"];
sz_prop_thresh = 0.01;

for a = 1:length(atlas)
    fprintf("Lausanne %d ", atlas(a))
    for s = 1:length(per_sz_or_all_sz)
        fprintf("results for %s ", per_sz_or_all_sz(s))
        for d = 1:length(det_method)
            fprintf("using %s onset detection \n", det_method(d))
            clear all_pat_table
            % Itearate across patients
            for pat = 1:size(final_output_both_atlas,1)
                patient = final_output_both_atlas.Patient_id{pat};
                pat_onset = final_output_both_atlas(string(final_output.Patient_id) == patient,: );
                if all(isnan(pat_onset.clo_chan{:})) | sum(pat_onset.clo_chan{:}) ==0
                    fprintf("%s does not have clinically labelled onset \n", patient)
                     continue
                end
                
                pat_clo = pat_onset.(sprintf("clo_roi_%d", atlas(a))){:};
                pat_auto = pat_onset.(sprintf("%s_roi_%d", det_method(d), atlas(a))){:};
                
                % Remove any seizures with onset in >50% of regions or no onset
                % detected
                rm_sz =  sum(pat_auto,1) == 0 |...
                    sum(pat_auto,1) >= size(pat_auto,1)/2;
                pat_auto = pat_auto(:,~rm_sz);
                if size(pat_auto,2) == 0
                    fprintf("%s all seizures removed \n", patient)
                    continue
                end
                pat_sz_ids = string(pat_onset.Segment_ids{:});
                pat_sz_id = pat_sz_ids(~rm_sz);
                pat_sz_types = pat_onset.sz_types{:}(~rm_sz);
                
                % Isolate the seizure types you're interesterested in
                if sz_types ~= "all"
                    keep_sz = ismember(pat_sz_types, sz_type);
                    pat_auto = pat_auto(:,keep_sz);
                    pat_sz_id = pat_onset.Segment_ids{:}(keep_sz);
                    % pat_sz_types = pat_onset.ilae_sz_type{:}(keep_sz);
                end
                
                % Compare each imprint onset against CLO (with permutation test)
                pat_comp_tab = comp_auto_clo(pat_clo,...
                    pat_auto, "per_sz_or_all_sz", per_sz_or_all_sz(s),...
                    "sz_prop_thresh", sz_prop_thresh);
                
                pat_surg_out = pat_onset.Surgery_outcome{:};
                pat_surg_year = pat_onset.("Surgery year");
                pat_out_year = pat_onset.("Outcome year"){:};
                year_1 = pat_out_year - pat_surg_year == 1;

                pat_comp_tab.patient_id = repmat(string(patient),size(pat_comp_tab,1),1);
                pat_comp_tab.outcome = repmat(pat_surg_out(year_1),size(pat_comp_tab,1),1);
                pat_comp_tab.op_type = repmat(pat_onset.("Op type"),size(pat_comp_tab,1),1);

                
                % Identify whether results are per seizure or summarised across seizures
                if per_sz_or_all_sz(s) == "per_sz"
                    pat_comp_tab.sz_id = pat_sz_id;
                        else
                    pat_comp_tab.sz_id = {pat_sz_id};
                end
        
                % Add patient to table
                if exist('all_pat_table', 'var')
                    all_pat_table = [all_pat_table; pat_comp_tab];
                else
                    all_pat_table = pat_comp_tab;
                end
                                   
                comp_table.(sprintf("laus_%d", atlas(a))).(sprintf("%s",per_sz_or_all_sz(s))).(sprintf("%s", det_method(d))) =...
                    all_pat_table;
            end
        end
    end
end
