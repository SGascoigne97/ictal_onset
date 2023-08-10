final_output = final_output(cellfun(@sum, final_output.resected_roi_120)~=0,:);
%%
bound_tbl = [repelem([0.01,0.1,0.2,0.25,1/3],7);...
    repmat([0.25,1/3,0.5,2/3,0.75,0.8,1],1,5)]';
chan_or_roi = "roi_120";

for bound = 3:size(bound_tbl,1)
    low = bound_tbl(bound,1);
    high = bound_tbl(bound,2);

    low_print = round(low*100);
    high_print = round(high*100);

    mkdir(sprintf('figures/rapid_prototype/comparing_consensus_thresh/brain_plots/low%d_high%d', low_print, high_print))
    
    for pat = 1:size(final_output,1)
        pat_onset = final_output(pat,:);
        imprint = pat_onset.imprint_roi_120{:};
        resec = double(pat_onset.resected_roi_120{:});
        
        prop_sz = sum(imprint,2)/size(imprint,2);
        in_high = double(prop_sz>=high);
        in_low = double(prop_sz>=low);
        in_low_not_high = in_low - in_high;
    
        pat_roi = pat_onset.roi_names_120{:};
    
        % See plotting limits to show relevant hemispheres
        if any(contains(pat_roi,"l.")) & any(contains(pat_roi,"r."))
            % Bilateral placement
            y_lim = [0 450];
            x_lim = [-500 1300];
    
        elseif any(contains(pat_roi,"l.")) & ~any(contains(pat_roi,"r."))
            % Left hemisphere placement
            y_lim = [230 450]; % If only left hemisphere
            x_lim = [-250 900];
    
        elseif ~any(contains(pat_roi,"l.")) & any(contains(pat_roi,"r."))
            % Right hemisphere placement
            y_lim = [0 210];
            x_lim = [-250 900];
        end
        
        if sum(in_low_not_high) ==0
            fprintf("%s no onset in %d<n<%d seizures \n", pat_onset.Patient_id{:}, low_print, high_print)
            continue
        end
    
        if sum(in_high) ==0
            fprintf("%s no onset in >=%d seizures \n", pat_onset.Patient_id{:}, high_print)
            continue
        end   
        
        fig = figure("Position", [0,0,1200,900], "Visible","off");
        subplot(4,2,1)
        cm = [.7,.7,.7; 1,0,0];
        if all(in_low == 1)
            continue
        end
        plotBrain_NE(pat_roi, in_low ,"cm", cm);
        title(sprintf("In >=%d%% of onsets", low_print))
        colorbar off
        ylim(y_lim)
        
        subplot(4,2,2)
        plotBrain_NE(pat_roi, in_high ,"cm", [.7,.7,.7; 1,0,0]);
        title(sprintf("In >=%d%% of onsets", high_print))
        colorbar off
        ylim(y_lim)
        
        subplot(4,1,2)
        plotBrain_NE(pat_roi, in_low_not_high ,"cm", [.7,.7,.7; 1,0,0]);
        clim([0,1])
        title(sprintf("In %d%%<=n<=%d%% of onsets", low_print, high_print))
        colorbar off
        ylim(y_lim)
        
        subplot(4,1,3)
        plotBrain_NE(pat_roi, resec ,"cm", [.7,.7,.7;1,1,0]);
        clim([0,1])
        title("Resected")
        colorbar off
        ylim(y_lim)
        xlim(x_lim)
        
        subplot(4,1,4)
        col_sch = (in_low_not_high*2)+resec; 
        % Resected only: yellow
        % In low not high, not resected: red
        % In low not high, resected: orange
        cm =  [.7,.7,.7;1,1,0;1,0,0;1,0.5,0];
        if max(col_sch) < size(cm,1)
            cm = cm(1:(max(col_sch)+1),:);    
        end
        cm = interp1(cm, 1:0.01:size(cm,1));
        plotBrain_NE(pat_roi, col_sch ,"cm", cm);
        colorbar off
        hold on
        cm_scatter = [1,1,0;1,0.5,0;1,0,0];
        for i= 1:3
            scatter(-100,-100, [], cm_scatter(i,:), 'filled')
        end
        hold off
        ylim(y_lim)
        xlim(x_lim)
        
        legend(["Resected only", "In some onsets, resec",...
            "In some onsets, not resec"], 'Location','northeast')
        
        legend boxoff  
        sgtitle(sprintf("%s, (ILAE %d)", pat_onset.Patient_id{:}, pat_onset.outcome))
        
        saveas(fig,sprintf('figures/rapid_prototype/comparing_consensus_thresh/brain_plots/low%d_high%d/%s.png', low_print, high_print, pat_onset.Patient_id{:}));
        
    end
end
