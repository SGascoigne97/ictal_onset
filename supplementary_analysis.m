%% Supplementary analyses

few_lim = 0.25;
most_lim = 0.75; 

comparison = "resection";
chan_or_roi = "roi";

mkdir('figures/across_patients/supplementary')
count2 = 1;
for det_method = ["imprint", "EI", "PLHG"]
    resec_roi_comp = readtable(sprintf('tables/across_pat_%s_%s_%s.csv', det_method, comparison, chan_or_roi));
    patients = unique(resec_roi_comp.Patient_id, 'stable');
    most_few_tab = array2table(zeros(length(unique(resec_roi_comp.Patient_id)),3));
    most_few_tab.Properties.VariableNames = {'Patient_ID', 'Prop_few', 'Prop_most'};
    most_few_tab.Patient_ID = patients;

    figure(count2)
    count = 1;
    for adj = -0.05:0.01:0.05
        for pat = 1:length(patients)
            pat_comp = resec_roi_comp(resec_roi_comp.Patient_id == patients(pat),:);
            most_few_tab.Prop_few(pat) = sum(pat_comp.Percentage_resec <=few_lim+adj)/size(pat_comp,1);
            most_few_tab.Prop_most(pat) = sum(pat_comp.Percentage_resec >=most_lim+adj)/size(pat_comp,1);
            most_few_tab.Outcome(pat) = discretize(pat_comp.Y1_outcome(1),[0 2.1 5.1], 'categorical',["Good", "Bad"]);
        end
    
        subplot(11,2,count)
        swarmchart( most_few_tab.Outcome, most_few_tab.Prop_few, 'filled')
        hold on
        boxchart( most_few_tab.Outcome, most_few_tab.Prop_few)
        hold off
        title(sprintf('Few (thresh = %.2f)', few_lim+adj))
        count = count + 1;
        subplot(11,2,count)
        swarmchart( most_few_tab.Outcome, most_few_tab.Prop_most, 'filled')
        hold on
        boxchart( most_few_tab.Outcome, most_few_tab.Prop_most)
        hold off
        title(sprintf('Most (thresh = %.2f)', most_lim+adj))
        count = count + 1;
       
    end
     sgtitle(sprintf("Testing thresholds (%s)", det_method))
     saveas(gcf,sprintf('figures/across_patients/supplementary/test_thresh_most_few_%s', det_method), 'png')
     count2 = count2+1;
end

