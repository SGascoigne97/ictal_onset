

sz_120 = comp_table.laus_120.per_sz.imprint;
sz_250 = comp_table.laus_250.per_sz.imprint;

for col_id = [1:4]
    sz_120.Properties.VariableNames(col_id) = ...
        sprintf("%s_120",string(sz_120.Properties.VariableNames(col_id)));
    sz_250.Properties.VariableNames(col_id) = ...
        sprintf("%s_250",string(sz_250.Properties.VariableNames(col_id)));
end

%join(sz_120, sz_250, 'Keys',"sz_type")
comb_tab = outerjoin(sz_120, sz_250, 'Type', 'Left', 'MergeKeys', true);

comb_tab.sz_type(string(comb_tab.sz_type) == 'f') = {'focal'};

figure()
subplot(221)
scatter(comb_tab.("Percentage of CLO captured_120"),...
    comb_tab.("Percentage of CLO captured_250"), 'filled')
title("Perc CLO")
xlabel("Laus 120")
ylabel("Laus 250")
lsline()
subplot(222)
gscatter(comb_tab.Perc_z_120, comb_tab.Perc_z_250,categorical(comb_tab.sz_type) ,'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Perc CLO (z)")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 8])
ylim([-2 8])
subplot(223)
scatter(comb_tab.Jaccard_120, comb_tab.Jaccard_250, 'filled')
lsline()
title("Jacc")
xlabel("Laus 120")
ylabel("Laus 250")
subplot(224)
gscatter(comb_tab.Jaccard_z_120, comb_tab.Jaccard_z_250,categorical(comb_tab.sz_type), 'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Jacc (z)")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 8])
ylim([-2 8])


%% 
comb_tab.op_type = categorical(comb_tab.op_type );
max_per_pat = groupsummary(comb_tab(:,[1:7,10:13]),["patient_id","op_type"],"max");

figure(col_id)
subplot(221)
gscatter(max_per_pat.max_Perc_z_120,max_per_pat.max_Perc_z_250 , categorical(max_per_pat.max_outcome>2),'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Perc CLO (z) by outcome (true = bad)")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 10])
ylim([-2 10])

subplot(222)
gscatter(max_per_pat.max_Jaccard_z_120, max_per_pat.max_Jaccard_z_250, categorical(max_per_pat.max_outcome>2), 'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Jacc (z) by outcome (true = bad)")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 10])
ylim([-2 10])

subplot(223)
gscatter(max_per_pat.max_Perc_z_120, max_per_pat.max_Perc_z_250, max_per_pat.op_type, 'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Perc CLO (z) by op type")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 10])
ylim([-2 10])

subplot(224)
gscatter(max_per_pat.max_Jaccard_z_120, max_per_pat.max_Jaccard_z_250, max_per_pat.op_type, 'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Jacc (z) by op type")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 10])
ylim([-2 10])


%%


max_per_pat = groupsummary(comb_tab(:,[1:7,9:13]),["patient_id","op_type","sz_type"],"max");
max_per_pat = max_per_pat(max_per_pat.sz_type ~= "unknown" &max_per_pat.sz_type ~= "sg" ,:);

figure()
subplot(211)
gscatter(max_per_pat.max_Perc_z_120, max_per_pat.max_Perc_z_250, max_per_pat.sz_type, 'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Perc CLO (z) by seizure type")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 10])
ylim([-2 10])

subplot(212)
gscatter(max_per_pat.max_Jaccard_z_120, max_per_pat.max_Jaccard_z_250, max_per_pat.sz_type, 'filled')
lsline()
xline(0, 'r')
yline(0, 'r')
title("Jacc (z) by seizure type")
xlabel("Laus 120")
ylabel("Laus 250")
xlim([-2 10])
ylim([-2 10])



%%
imp_120 = comp_table.laus_120.per_sz.imprint;
ei_120 = comp_table.laus_120.per_sz.EI;
plhg_120 =  comp_table.laus_120.per_sz.PLHG;

for col_id = [1:4]
    imp_120.Properties.VariableNames(col_id) = ...
        sprintf("%s_imp",string(imp_120.Properties.VariableNames(col_id)));
    ei_120.Properties.VariableNames(col_id) = ...
        sprintf("%s_ei",string(ei_120.Properties.VariableNames(col_id)));
     plhg_120.Properties.VariableNames(col_id) = ...
        sprintf("%s_plhg",string(plhg_120.Properties.VariableNames(col_id)));
end


imp_vs_ei = outerjoin(imp_120, ei_120, 'Type', 'Left', 'MergeKeys', true);
imp_ei_plhg = outerjoin(imp_vs_ei, plhg_120, 'Type', 'Left', 'MergeKeys', true);

figure()
subplot(221)
histogram(imp_ei_plhg.Perc_z_imp-imp_ei_plhg.Perc_z_ei)
xline(0, 'r')
title("Imprint - EI % CLO (z)")
ylim([0 100])

subplot(222)
histogram(imp_ei_plhg.Perc_z_imp-imp_ei_plhg.Perc_z_plhg)
xline(0, 'r')
title("Imprint - PLHG % CLO (z)")
ylim([0 100])

subplot(223)
histogram(imp_ei_plhg.Jaccard_z_imp-imp_ei_plhg.Jaccard_z_ei)
xline(0, 'r')
title("Imprint - EI Jacc (z)")
ylim([0 100])

subplot(224)
histogram(imp_ei_plhg.Jaccard_z_imp-imp_ei_plhg.Jaccard_z_plhg)
xline(0, 'r')
title("Imprint - PLHG Jacc (z)")
ylim([0 100])

%% Plot by seizure type
subclin = imp_ei_plhg(imp_ei_plhg.sz_type == "subclin",:);
focal = imp_ei_plhg(imp_ei_plhg.sz_type == "focal",:);

figure()
subplot(221)
histogram(subclin.Perc_z_imp-subclin.Perc_z_ei, BinWidth=0.5, FaceColor='g' )
hold on
histogram(focal.Perc_z_imp-focal.Perc_z_ei, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
title("Imprint - EI % CLO (z)")
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

subplot(2,2,2)
histogram(subclin.Jaccard_z_imp-subclin.Jaccard_z_ei, BinWidth=0.5, FaceColor='g' )
hold on
histogram(focal.Jaccard_z_imp-focal.Jaccard_z_ei, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
title("Imprint - EI Jacc (z)")
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

subplot(2,2,3)
histogram(subclin.Perc_z_imp-subclin.Perc_z_plhg, BinWidth=0.5, FaceColor='g' )
hold on
histogram(focal.Perc_z_imp-focal.Perc_z_plhg, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
title("Imprint - PLHG % CLO (z)")
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

subplot(2,2,4)
histogram(subclin.Jaccard_z_imp-subclin.Jaccard_z_plhg, BinWidth=0.5, FaceColor='g' )
hold on
histogram(focal.Jaccard_z_imp-focal.Jaccard_z_plhg, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
title("Imprint - PLHG Jacc (z)")
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})




%%
max_per_pat = groupsummary(imp_ei_plhg(:,[1:6,9:end]),["patient_id","sz_type"],"max");
subclin = max_per_pat(max_per_pat.sz_type == "subclin",:);
focal = max_per_pat(max_per_pat.sz_type == "focal",:);
figure()
subplot(2,2,1)
histogram(subclin.max_Perc_z_imp-subclin.max_Perc_z_ei, BinWidth=0.5, FaceColor='g')
title("Imprint - EI % CLO (z)")
hold on
histogram(focal.max_Perc_z_imp-focal.max_Perc_z_ei, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

subplot(2,2,2)
histogram(subclin.max_Jaccard_z_imp-subclin.max_Jaccard_z_ei, BinWidth=0.5, FaceColor='g')
title("Imprint - EI Jacc (z)")
hold on
histogram(focal.max_Jaccard_z_imp-focal.max_Jaccard_z_ei, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

subplot(2,2,3)
histogram(subclin.max_Perc_z_imp-subclin.max_Perc_z_plhg, BinWidth=0.5, FaceColor='g')
title("Imprint - PLHG % CLO (z)")
hold on
histogram(focal.max_Perc_z_imp-focal.max_Perc_z_plhg, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

subplot(2,2,4)
histogram(subclin.max_Jaccard_z_imp-subclin.max_Jaccard_z_plhg, BinWidth=0.5, FaceColor='g')
title("Imprint - PLHG Jacc (z)")
hold on
histogram(focal.max_Jaccard_z_imp-focal.max_Jaccard_z_plhg, BinWidth=0.5, FaceAlpha=0.5, FaceColor='r')
hold off
xline(0, LineWidth=2)
xlim([-6 6])
legend({"Subclin", "Focal"})

%%
figure()
subplot(2,2,1)
histogram(max_per_pat.max_Perc_z_imp-max_per_pat.max_Perc_z_ei, BinWidth=0.5)
title("Imprint - EI % CLO (z)")
xline(0, LineWidth=2)
xlim([-6 6])

subplot(2,2,2)
histogram(max_per_pat.max_Jaccard_z_imp-max_per_pat.max_Jaccard_z_ei, BinWidth=0.5)
title("Imprint - EI Jacc (z)")
xline(0, LineWidth=2)
xlim([-6 6])

subplot(2,2,3)
histogram(max_per_pat.max_Perc_z_imp-max_per_pat.max_Perc_z_plhg, BinWidth=0.5)
title("Imprint - PLHG % CLO (z)")
xline(0, LineWidth=2)
xlim([-6 6])

subplot(2,2,4)
histogram(max_per_pat.max_Jaccard_z_imp-max_per_pat.max_Jaccard_z_plhg, BinWidth=0.5)
title("Imprint - PLHG Jacc (z)")
xline(0, LineWidth=2)
xlim([-6 6])

