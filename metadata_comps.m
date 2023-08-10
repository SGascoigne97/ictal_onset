%% Look into metadata (to report in supplementary)

[tab_out, ~, ~, lab_out] = crosstab(final_output_all.outcome);
[tab_op, ~, ~, lab_op] = crosstab(final_output_all.("Op type"));

surg = string(final_output_all.("Op type"));
surg(surg ~= "F Lx" & surg ~= "T Lx") = "other";
figure()
boxchart(categorical(surg),final_output_all.outcome)
hold on
swarmchart(categorical(surg),final_output_all.outcome, 'filled')
hold off
ylim([0,6])
ylabel("ILAE outcome")
xlabel("Surgery type")

%%
final_output_sz_types = final_output_all(cellfun(@iscell, final_output_all.sz_types),:);
sz_types = cellfun(@unique, final_output_sz_types.sz_types, 'UniformOutput', false);

figure(2)
swarmchart(cellfun(@length, sz_types), final_output_sz_types.outcome, 'filled')
hold on
boxchart(cellfun(@length, sz_types), final_output_sz_types.outcome)
hold off
ylim([0,6])

%% Separate into focal only, subclinical only, and mixed

%% Separate into sg present or not
for pat = 1:size(final_output_all)
    sz_types = string(final_output_all(pat,:).sz_types{:});
    has_ftbtc(pat) = contains("sg", sz_types);
end
figure(3)
swarmchart(double(has_ftbtc), final_output_all.outcome, 'filled')
hold on
boxchart(double(has_ftbtc), final_output_all.outcome)
hold off
ylim([0,6])
ylabel("ILAE outcome")
xlabel("Presence of FTBTC")
set(gca, "XTick", [0,1], "XTickLabel", ["No", "Yes"])

figure(1)
subplot(2,2,1)
boxchart(categorical(final_output.sex), final_output.outcome)
hold on
swarmchart(categorical(final_output.sex), final_output.outcome, 'filled')
hold off
xlabel("Sex")
ylabel("Outcome")

subplot(2,2,2)
boxchart(double(final_output.outcome>3), final_output.age, "MarkerStyle","none")
hold on
swarmchart(double(final_output.outcome>3),final_output.age, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
xlabel("Outcome")
ylabel("Age (years)")

subplot(2,2,3)
boxchart(double(final_output.outcome>3), (final_output.age-final_output.onset_age), "MarkerStyle","none")
hold on
swarmchart(double(final_output.outcome>3),(final_output.age-final_output.onset_age), 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
xlabel("Outcome")
ylabel("Epilepsy Duration (years)")

subplot(2,2,4)
gscatter(final_output.age, (final_output.age-final_output.onset_age), final_output.outcome>3, 'filled' )
xlabel("Age (years)")
ylabel("Epilepsy Duration (years)")
lsline()
legend(["ILAE 1-2", "ILAE 3+"], "Location", "northwest")
%%
final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
prop_resec = nan(size(final_output_all,1),1);
for pat = 1:size(final_output_all,1)
    resec = final_output_all(pat,:).resected_roi_120{:};
    prop_resec(pat) = sum(resec)/length(resec);
end


out = final_output_all.outcome;
% Remove patients with no recorded regions resected;
out = out(prop_resec>0);
prop_resec = prop_resec(prop_resec>0);


figure(1)
swarmchart(double(out>2), prop_resec,'filled')
hold on
boxchart(double(out>2), prop_resec)
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
ylabel("Proportion of recorded regions resected")

%%
prop_clo = nan(size(final_output_all,1),1);
for pat = 1:size(final_output_all,1)
    clo = final_output_all(pat,:).clo_roi_120{:};
    prop_clo(pat) = sum(clo)/length(clo);
end


out = final_output_all.outcome;
% Remove patients with no recorded regions resected;
out = out(prop_clo>0);
prop_clo = prop_clo(prop_clo>0);

figure(2)
swarmchart(double(out>2), prop_clo,'filled')
hold on
boxchart(double(out>2), prop_clo)
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
ylabel("Proportion of recorded regions labelled as CLO")

%%
figure(3)
tiledlayout(1,3)
for chan_or_roi = ["chan" "roi_120", "roi_250"]
    nexttile
    prop_resec = nan(size(final_output_all,1),1);
    for pat = 1:size(final_output_all,1)
        resec = final_output_all(pat,:).(sprintf("resected_%s", chan_or_roi)){:};
        prop_resec(pat) = sum(resec)/length(resec);
    end
    prop_clo = nan(size(final_output_all,1),1);
    for pat = 1:size(final_output_all,1)
        clo = final_output_all(pat,:).(sprintf("clo_%s", chan_or_roi)){:};
        prop_clo(pat) = sum(clo)/length(clo);
    end
    out = final_output_all.outcome;
    
    gscatter(prop_resec, prop_clo, out<3, 'filled')
    lsline()
    xlim([0,1])
    ylim([0,1])
    axis square
    xlabel("Proportion resected")
    ylabel("Proportion clo")
    legend({"ILAE 1-2", "ILAE 3+"})
end

