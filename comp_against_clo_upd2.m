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

atlas = [120, 250];
per_sz_or_all_sz = ["all_sz", "per_sz"];
det_method = ["imprint", "EI", "PLHG"];

% Paramaters to select at start (will not be scanned)
sz_types = "all"; %["focal", "subclin", "all"];
sz_prop_thresh = 0.01;


% output_tab = load('tables/final_output.mat');
% final_output = output_tab.final_output;
% clear output_tab

% Remove those without CLO
rm_pat = cellfun(@sum,final_output.clo_chan) == 0;
final_output = final_output(~rm_pat,:);
%%

for a = 1:length(atlas)
    fprintf("Lausanne %d ", atlas(a))
    for s = 1%:length(per_sz_or_all_sz)
        fprintf("results for %s ", per_sz_or_all_sz(s))
        for d = 1:length(det_method)
            fprintf("using %s onset detection \n", det_method(d))
            clear all_pat_table
            % Itearate across patients
            for pat = 1:size(final_output,1)
                patient = final_output.Patient_id{pat};
                pat_onset = final_output(string(final_output.Patient_id) == patient,: );
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
                    pat_sz_types = pat_onset.ilae_sz_type{:}(keep_sz);
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
                    pat_comp_tab.sz_type = pat_sz_types;
                else
                    pat_comp_tab.sz_id = {pat_sz_id};
                    pat_comp_tab.sz_type = {pat_sz_types};
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

%%
param_tab = table(repmat(repelem([120;250],2,1),3,1),...
    repelem(["imprint";"EI";"PLHG"],4,1), repmat(["per_sz";"all_sz"],6,1),...
    'VariableNames', ["atl", "det", "per_or_all"]);
sz_prop_thresh = 0.01;

% param_tab = param_tab(2:2:12,:);
param_tab = param_tab([1,3],:);
for r = 1:size(param_tab,1)
    param_output = onset_comp_clo(final_output, "atl",param_tab.atl(r),...
        "per_or_all",param_tab.per_or_all(r),"det",param_tab.det(r),...
        "sz_prop_thresh", sz_prop_thresh, "sz_types", "all");
     comp_table.(sprintf("laus_%d", param_tab.atl(r))).(sprintf("%s",param_tab.per_or_all(r))).(sprintf("%s", param_tab.det(r))) =...
            param_output;
end

%% 
columns = ["patient_id", "sz_id", "Perc_z", "Jaccard_z", "outcome", "op_type", "sz_type"];
imprint_tab = comp_table.laus_120.per_sz.imprint(:,columns);
ei_tab = comp_table.laus_120.per_sz.EI(:,columns);
plhg_tab = comp_table.laus_120.per_sz.PLHG(:,columns);

for col = 3:4
    imprint_tab.Properties.VariableNames(col) = sprintf("%s_imprint", string(imprint_tab.Properties.VariableNames(col)));
    ei_tab.Properties.VariableNames(col) = sprintf("%s_ei", string(ei_tab.Properties.VariableNames(col)));
    plhg_tab.Properties.VariableNames(col) = sprintf("%s_plhg", string(plhg_tab.Properties.VariableNames(col)));
end

%join(imprint_tab, ei_tab, "Keys", ["patient_id", "sz_id", "sz_type", "op_type", "outcome"])
final_tab = outerjoin(imprint_tab, ei_tab, 'Type', 'Left', 'MergeKeys', true);
final_tab = outerjoin(final_tab, plhg_tab, 'Type', 'Left', 'MergeKeys', true);

final_tab = final_tab(~isnan(final_tab.Perc_z_ei),:);
%%
figure(1)
tiledlayout(3,2)
method_mat = ["imprint", "ei"; "imprint", "plhg"; "ei", "plhg"];
for combinat = 1:3
    nexttile
    gscatter(final_tab.(sprintf("Perc_z_%s", method_mat(combinat,1))),...
        final_tab.(sprintf("Perc_z_%s", method_mat(combinat,2))), ...
        [], 'filled')
    lsline()
    title(sprintf("Percentage of CLO z %s vs %s",method_mat(combinat,1),method_mat(combinat,2)))
    xlabel(method_mat(combinat,1))
    ylabel(method_mat(combinat,2))
    xline(0, 'k')
    yline(0, 'k')
    xlim([-2 8])
    ylim([-2 8])
    nexttile
    gscatter(final_tab.(sprintf("Jaccard_z_%s", method_mat(combinat,1))),...
        final_tab.(sprintf("Jaccard_z_%s", method_mat(combinat,2))), ...
        [], 'filled')
    lsline()
    title(sprintf("Jaccard z %s vs %s",method_mat(combinat,1),method_mat(combinat,2)))
    xlabel(method_mat(combinat,1))
    ylabel(method_mat(combinat,2))
    xline(0, 'k')
    yline(0, 'k')
    xlim([-2 8])
    ylim([-2 8])
end

%%
figure(2)
tiledlayout(3,2)
method_mat = ["imprint", "ei"; "imprint", "plhg"; "ei", "plhg"];
for combinat = 1:3
    nexttile
    histogram(final_tab.(sprintf("Perc_z_%s", method_mat(combinat,1)))...
        -final_tab.(sprintf("Perc_z_%s", method_mat(combinat,2))), BinWidth=0.5)
    title(sprintf("Percentage of CLO z %s - %s",method_mat(combinat,1),method_mat(combinat,2)))
    xline(-1.96, 'k')
    xline(1.96, 'k')
    nexttile
    histogram(final_tab.(sprintf("Jaccard_z_%s", method_mat(combinat,1)))-...
        final_tab.(sprintf("Jaccard_z_%s", method_mat(combinat,2))), BinWidth=0.5)
    title(sprintf("Jaccard z %s - %s",method_mat(combinat,1),method_mat(combinat,2)))
 
    xline(-1.96, 'k')
    xline(1.96, 'k')
end

%%
figure()
subplot(131)
imagesc(final_output.imprint_chan{1,1})
title("Channels")
set(gca, "YTick", 1:length(final_output.channel_names{1,1}), "YTickLabel",final_output.channel_names{1,1} )
subplot(132)
imagesc(final_output.imprint_roi_120{1,1})
title("Lausanne 120")
set(gca, "YTick", 1:length(final_output.roi_names_120{1,1}), "YTickLabel",final_output.roi_names_120{1,1} )
subplot(133)
imagesc(final_output.imprint_roi_250{1,1})
title("lausanne 250")
set(gca, "YTick", 1:length(final_output.roi_names_250{1,1}), "YTickLabel",final_output.roi_names_250{1,1} )
sgtitle(final_output.Patient_id(1))

%%
tab_120 = comp_table.laus_120.per_sz.imprint;
tab_250 = comp_table.laus_250.per_sz.imprint;

robust_imp_tab = innerjoin(tab_120, tab_250, Keys="sz_id");

robust_imp_tab.sz_type = robust_imp_tab.sz_type_tab_120;
robust_imp_tab(robust_imp_tab.sz_type_tab_120 == "f",:).sz_type=...
    repelem({"focal"},length(robust_imp_tab(robust_imp_tab.sz_type_tab_120 == "f",:).sz_type_tab_120),1);
robust_imp_tab(string(robust_imp_tab.sz_type) == "u",:).sz_type =...
    repelem({"unknown"},length(robust_imp_tab(string(robust_imp_tab.sz_type) == "u",:).sz_type),1);


%%
figure(1)
subplot(121)
gscatter(robust_imp_tab.Perc_z_tab_120, robust_imp_tab.Perc_z_tab_250, string(robust_imp_tab.sz_type) , 'filled')
set(gca, "XLim", [-3,8], "YLim", [-3,8])
xlabel("Lausanne 120")
ylabel("Lausanne 250")
lsline()
xline(0, '--')
yline(0, '--')
title("Percentage of CLO (z)")
subplot(122)
gscatter(robust_imp_tab.Jaccard_z_tab_120, robust_imp_tab.Jaccard_z_tab_250,string(robust_imp_tab.sz_type) , 'filled')
set(gca, "XLim", [-3,8], "YLim", [-3,8])
xlabel("Lausanne 120")
ylabel("Lausanne 250")
lsline()
xline(0, '--')
yline(0, '--')
title("Jaccard (z)")


figure(2)
subplot(221)
gscatter(robust_imp_tab.Perc_z_tab_120, robust_imp_tab.Perc_z_tab_250, string(robust_imp_tab.op_type_tab_250) , 'filled')
set(gca, "XLim", [-3,8], "YLim", [-3,8])
xlabel("Lausanne 120")
ylabel("Lausanne 250")
lsline()
xline(0, '--')
yline(0, '--')
title("Percentage of CLO (z)")
subplot(222)
gscatter(robust_imp_tab.Jaccard_z_tab_120, robust_imp_tab.Jaccard_z_tab_250,string(robust_imp_tab.op_type_tab_250) , 'filled')
set(gca, "XLim", [-3,8], "YLim", [-3,8])
xlabel("Lausanne 120")
ylabel("Lausanne 250")
lsline()
xline(0, '--')
yline(0, '--')
title("Jaccard (z)")

subplot(223)
gscatter(robust_imp_tab.Perc_z_tab_120, robust_imp_tab.Perc_z_tab_250, string(robust_imp_tab.outcome_tab_120) , 'filled')
set(gca, "XLim", [-3,8], "YLim", [-3,8])
xlabel("Lausanne 120")
ylabel("Lausanne 250")
lsline()
xline(0, '--')
yline(0, '--')
title("Percentage of CLO (z)")
subplot(224)
gscatter(robust_imp_tab.Jaccard_z_tab_120, robust_imp_tab.Jaccard_z_tab_250,string(robust_imp_tab.outcome_tab_120) , 'filled')
set(gca, "XLim", [-3,8], "YLim", [-3,8])
xlabel("Lausanne 120")
ylabel("Lausanne 250")
lsline()
xline(0, '--')
yline(0, '--')
title("Jaccard (z)")
