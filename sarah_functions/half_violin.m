function [] = half_violin(main_data_tab, det_meths, comp_measures, opts)
% Compute average EEG across regions of interest for each patient

% input:
%   - data_tbl: full data table
%   - optional inputs
%       - sz_type: seizure type of interest
%       - min_sz_duration: minimum duration for each seizure 
%       - min_sz_count: minimum number of seizures to ahve been recorded per patient

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        main_data_tab
        det_meths
        comp_measures
        opts.parc = "roi_120";
        opts.onset_across = 1;
        opts.summ_meas = [];
        opts.vis_plot = "on"
        opts.save_fig = 0;
        opts.save_loc = 'figures/'
        opts.file_type = "png"

    end
    
    %fill in optional arguments
    parc = opts.parc;
    onset_across = opts.onset_across;
    summ_meas = opts.summ_meas;
    vis_plot = opts.vis_plot;
    save_fig = opts.save_fig;
    save_loc = opts.save_loc;
    file_type = opts.file_type;

    
    onset_across_titles = [sprintf("Across subject %s",summ_meas), "Onset across seizures"];
    % Create figure 
    fig = figure("Position", [10, 10, 1800, 900], 'Visible',vis_plot);
    
    if length(comp_measures)>1
        tiledlayout(ceil(length(comp_measures))/2,2)
    else
        tiledlayout(1,1)
    end

    
    rng(6)
    for comp_measure = comp_measures
        if comp_measure == "Perc" | comp_measure == "Jacc"
            y_lim = [-0.1,1.1];
            %bin_wd = 0.1;
        elseif comp_measure == "Coh"
            y_lim = [-0.6,1.1];
            %bin_wd = 0.1;
        elseif contains(comp_measure, "z")
            y_lim = [-3, 12]; % May need to tweak this later
            %bin_wd = 1;
        end
    
        % Compute offset between violin plots as the maximum number of subjects
        % across methods
        subjs = [];
        for det_meth = det_meths
            subjs = [subjs, size(main_data_tab.across_sz.(sprintf(parc)).(sprintf(det_meth)),1)];
        end
        offset = (max(subjs)+10)/2;
        max_x = (offset*((length(det_meths))-1)) +10;
        x_lim = [-10, max_x];
        
        nexttile
        hold on
        for det_ind = 1:length(det_meths)
            if onset_across == 1
                data_tab =  main_data_tab.across_sz.(sprintf("%s", parc)).(sprintf("%s", det_meths(det_ind)));
                data = data_tab.(sprintf("%s", comp_measure));
                out = data_tab.Outcome;
            else
                data_tab = main_data_tab.per_sz.(sprintf("%s", parc)).(sprintf("%s", det_meths(det_ind))).(sprintf("%s", summ_meas));
                data = data_tab.(sprintf("%s", comp_measure));
                out = data_tab.Outcome;
            end
            
            out = out(~isnan(data));
            data = data(~isnan(data));

            bin_wd = 5*(max(data)-min(data))/length(data);
           
            g = data(out<3);
            if isempty(g)
                continue
            end
            b = data(out>2);
            if isempty(b)
                continue
            end
            
            min_val_g = round(min(g)*(1/bin_wd))*bin_wd;
            min_val_b = round(min(b)*(1/bin_wd))*bin_wd;
            max_val_g = round(max(g)*(1/bin_wd))*bin_wd;
            max_val_b = round(max(b)*(1/bin_wd))*bin_wd;

            if min_val_g > min(g)
                min_val_g = min_val_g - bin_wd;
            end 
            if min_val_b > min(b)
                min_val_b = min_val_b - bin_wd;
            end
            if max_val_g < max(g)
                max_val_g = max_val_g + bin_wd;
            end
            if max_val_b < max(b)
                max_val_b = max_val_b + bin_wd;
            end

            if min_val_g == max_val_g
                max_val_g = max_val_g+bin_wd;
            end
            if min_val_b == max_val_b
                max_val_b = max_val_b+bin_wd;
            end
       
            hist_vals_g = histcounts(g, min_val_g:bin_wd:max_val_g);
            % Add scatter points
            non_zero_grps = hist_vals_g(hist_vals_g~=0);
            jitt = [];
            for grp_sz = non_zero_grps
                jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
            end
            scatter((det_ind-1)*offset+jitt, sort(g), [], [0.3010 0.7450 0.9330], 'filled');
        
            hist_vals_b = histcounts(b, min_val_b:bin_wd:max_val_b);
            
            % Add scatter points
            non_zero_grps = hist_vals_b(hist_vals_b~=0);
            jitt = [];
            for grp_sz = non_zero_grps
                jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
            end
            scatter((det_ind-1)*offset-jitt, sort(b), [],[0.8500 0.3250 0.0980], 'filled');
            if length(hist_vals_g) >1
                smooth_data_g = smoothdata(interp1(hist_vals_g, 1:0.1:length(hist_vals_g)));
                fill((det_ind-1)*offset+[0 smooth_data_g 0], ...
                    rescale((1:length([0 smooth_data_g 0]))/length([0 smooth_data_g 0]),...
                    min(g) - (bin_wd/10), max(g) + (bin_wd/10)), [0.3010 0.7450 0.9330], 'FaceAlpha',0.3,...
                    'LineStyle','none')
            end 

            if length(hist_vals_b) > 1
                smooth_data_b = smoothdata(interp1(hist_vals_b, 1:0.1:length(hist_vals_b)));
                fill((det_ind-1)*offset -[0 smooth_data_b 0],...
                    rescale((1:length([0 smooth_data_b 0]))/length([0 smooth_data_b 0]),...
                    min(b) - (bin_wd/10), max(b) + (bin_wd/10)), [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,...
                    'LineStyle','none')
            end
        
        end
        set(gca, "XTick", (0:(length(det_meths))-1)*offset, "XTickLabel", det_meths)
        title(comp_measure)
        ylim(y_lim)
        xlim(x_lim)
    
        if contains(comp_measure, "z")
            patch([x_lim(1), x_lim(2), x_lim(2), x_lim(1)], [-2, -2, 2, 2], [0.7,0.7,0.7], "FaceAlpha", 0.3, "LineStyle", "none")
        elseif comp_measure == "Coh"
            yline(0.2)
        end
    end
    sgtitle(sprintf("%s-wise comparison against resection (%s) ", parc, onset_across_titles(onset_across+1)))

    if save_fig == 1
        saveas(fig,sprintf("%sviolin_%s_%s.%s", save_loc, parc, strrep(onset_across_titles(onset_across+1), ' ', '_'),file_type))
    end
end
  
