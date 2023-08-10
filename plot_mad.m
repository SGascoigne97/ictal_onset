pat_dir = dir('MAD_scan_tabs');

%%
clear all_pat_mad
for pat = 3:size(pat_dir,1)
    pat_id = string(pat_dir(pat,:).name); % Is an issue when partient IDs arent the same length 
    load(sprintf('MAD_scan_tabs/%s', pat_id))
    if size(pat_tab,2) == 302
        if exist("all_pat_mad", "var")
            all_pat_mad = [all_pat_mad; pat_tab];
        else 
            all_pat_mad = pat_tab;
        end
    elseif size(pat_tab,2) == 102
         if exist("second_all_pat_mad", "var")
            second_all_pat_mad = [second_all_pat_mad; pat_tab];
        else 
            second_all_pat_mad = pat_tab;
         end
    else
        sprintf("%s not enough MAD values", pat_id)

    end
end
%%
%all_pat_mad_third = all_pat_mad(randi(size(all_pat_mad,1), 1,size(all_pat_mad,1)/3),:);
% 

tau = 5;
comp_vars =["prop_act", "prop_chan_ons", "ons_time"];
max_x = 15;

for comp = 1:length(comp_vars)
    comp_var = comp_vars(comp);

    mad_cols = contains(all_pat_mad.Properties.VariableNames,sprintf("%s_mad_", comp_var));
    mad_val = str2double(extractAfter(string(all_pat_mad.Properties.VariableNames(mad_cols)),"mad_"));
    mad_mat =table2array(all_pat_mad(:,mad_cols));
    diff_mat = mad_mat(:,1:end-1) - mad_mat(:,2:end);
    
    % Format plot titles and axis limits
    if comp_var == "prop_act"
        titles = {["A. Proportion of windows with activity detected", "across \tau values for all seizures"],...
            ["B. Change in proportion of windows with activity", "between \tau_i and \tau_{i+1} for all seizures"],...
            ["C. Mean and SD of proportion of windows ", "with activity across \tau values"],...
            ["D. Mean and SD of change in proportion of windows", " with activity between \tau_i and \tau_{i+1}"]};
        ylim_raw = [0,1];
        ylim_diff = [-0.2, 0.6];
    elseif comp_var == "prop_chan_ons"
         titles = {["A. Proportion of channels in onset", "across \tau values for all seizures"],...
            ["B. Change in proportion of channels in onset", "between \tau_i and \tau_{i+1} for all seizures"],...
            ["C. Mean and SD of proportion of channels ", "in onset across \tau values"],...
            ["D. Mean and SD of change in proportion of channels", " in onset between \tau_i and \tau_{i+1}"]};
         ylim_raw = [0,1];
         ylim_diff = [-0.6, 0.6];
    else 
          titles = {["A. Onset time (in windows)", "across \tau values for all seizures"],...
            ["B. Change in onset time between \tau_i and", "\tau_{i+1} for all seizures"],...
            ["C. Mean and SD of onset time (in windows)" , "across \tau values"],...
            ["D. Mean and SD of change in onset time", " between \tau_i and \tau_{i+1}"]};
         ylim_raw = [0,1500];
         ylim_diff = [-800,10];
    end

    xlim_all = [0, max_x];
    
    figure(comp)
    subplot(221)
    plot(mad_val, mad_mat);
    title(titles{1})
    xline(tau, '--r')
    yline(0, '--k')
    ylim(ylim_raw)
    xlim(xlim_all)
    
    subplot(222)
    plot(mad_val(2:end), diff_mat)
    title(titles{2})
    xline(tau, '--r')
    ylim(ylim_diff)
    xlim(xlim_all)
    
    mean_mad = mean(mad_mat, 'omitnan');
    std_mad = std(mad_mat, 'omitnan');
    subplot(223)
    plot(mad_val, mean_mad,'b', 'LineWidth', 2);
    hold on
    plot(mad_val, mean_mad + std_mad, ':b');
    plot(mad_val, mean_mad - std_mad, ':b');
    hold off
    title(titles{3})
    xline(tau, '--r')
    yline(0, '--k')
    ylim(ylim_raw)
    xlim(xlim_all)

    if comp_var == "ons_time"
        ylim_diff = [-100,100];
    end
    
    subplot(224)
    mean_diff = mean(diff_mat, 'omitnan');
    std_diff = std(diff_mat, 'omitnan');
    plot(mad_val(2:end), mean_diff, 'b', 'LineWidth', 2);
    hold on
    plot(mad_val(2:end), mean_diff + std_diff, ':b');
    plot(mad_val(2:end), mean_diff - std_diff, ':b');
    hold off
    title(titles{4})
    xline(tau, '--r')
    ylim(ylim_diff)
    xlim(xlim_all)

end



% mad_cols = contains(all_pat_mad.Properties.VariableNames,"mad_");
% mad_val = str2double(extractAfter(string(all_pat_mad.Properties.VariableNames(mad_cols)),"mad_"));
% mad_mat =table2array(all_pat_mad(:,mad_cols));
% diff_mat = mad_mat(:,1:end-1) - mad_mat(:,2:end);
% 
% figure(1)
% subplot(221)
% plot(mad_val, mad_mat);
% title(["A. Proportion of windows with activity detected", "across \tau values for all seizures"])
% xline(4, '--r')
% 
% subplot(222)
% plot(mad_val(2:end), diff_mat)
% title(["B. Change in proportion of windows with activity", "between \tau_i and \tau_{i+1} for all seizures"])
% xline(4, '--r')
% 
% subplot(223)
% plot(mad_val, mean(mad_mat),'b', 'LineWidth', 2);
% hold on
% plot(mad_val, mean(mad_mat) + std(mad_mat), ':b');
% plot(mad_val, mean(mad_mat) - std(mad_mat), ':b');
% hold off
% title(["C. Mean and SD of proportion of windows ", "with activity across \tau values"])
% 
% xline(4, '--r')
% 
% subplot(224)
% mean_diff = mean(diff_mat);
% std_diff = std(diff_mat);
% plot(mad_val(2:end), mean_diff, 'b', 'LineWidth', 2);
% hold on
% plot(mad_val(2:end), mean_diff + std_diff, ':b');
% plot(mad_val(2:end), mean_diff - std_diff, ':b');
% hold off
% title(["D. Mean and SD of change in proportion of windows", " with activity between \tau_i and \tau_{i+1}"])
% xline(4, '--r')