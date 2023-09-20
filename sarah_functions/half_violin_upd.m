function [] = half_violin_upd(data_tab, comp_var, grps, opts)
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
        data_tab table
        comp_var string
        grps (1,:) double
        opts.y_lim (1,2) double = [NaN,NaN];
        opts.grp_names = ["Group One", "Group Two"];
        opts.save_fig = 0;
        opts.save_loc = "figures/"
        opts.file_type = "png"

    end
    
    %fill in optional arguments
    y_lim = opts.y_lim;
    save_fig = opts.save_fig;
    save_loc = opts.save_loc;
    file_type = opts.file_type;
    grp_names = opts.grp_names;

    grp_vals = unique(grps);

    offset = (size(data_tab,1)+10)/2;
    max_x = offset + 10;
    x_lim = [-10, max_x];
    
    f = figure();
    f.Position = [500,500,200*length(comp_var),500];
    tiledlayout(1,length(comp_var)+1)
    
    var_ind = 1;
    for var = comp_var
        data_tab_var = data_tab(~isnan(data_tab.(sprintf(var))),:);
        grps_var = grps(~isnan(data_tab.(sprintf(var))));
        nexttile
        data = data_tab_var.(sprintf(var));

        hold on
    
        bin_wd = 5*range(data)/length(data);
        grp_a = data(grps_var == grp_vals(1));
        grp_b = data(grps_var == grp_vals(2));
                  
        min_val_a = round(min(grp_a)*(1/bin_wd))*bin_wd;
        min_val_b = round(min(grp_b)*(1/bin_wd))*bin_wd;
        max_val_a = round(max(grp_a)*(1/bin_wd))*bin_wd;
        max_val_b = round(max(grp_b)*(1/bin_wd))*bin_wd;
        
        if min_val_a > min(grp_a)
            min_val_a = min_val_a - bin_wd;
        end 
        if min_val_b > min(grp_b)
            min_val_b = min_val_b - bin_wd;
        end
        if max_val_a < max(grp_a)
            max_val_a = max_val_a + bin_wd;
        end
        if max_val_b < max(grp_b)
            max_val_b = max_val_b + bin_wd;
        end
        
        if min_val_a == max_val_a
            max_val_a = max_val_a+bin_wd;
        end
        if min_val_b == max_val_b
            max_val_b = max_val_b+bin_wd;
        end
        
        hist_vals_a = histcounts(grp_a, min_val_a:bin_wd:max_val_a);
        % Add scatter points
        non_zero_grps = hist_vals_a(hist_vals_a~=0);
        jitt = [];
        for grp_sz = non_zero_grps
            jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
        end
        scatter((var_ind-1)*offset+jitt, sort(grp_a), [], [0.3010 0.7450 0.9330], 'filled');
        
        hist_vals_b = histcounts(grp_b, min_val_b:bin_wd:max_val_b);
        % Add scatter points
        non_zero_grps = hist_vals_b(hist_vals_b~=0);
        jitt = [];
        for grp_sz = non_zero_grps
            jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
        end
        scatter((var_ind-1)*offset-jitt, sort(grp_b), [],[0.8500 0.3250 0.0980], 'filled');
        
        if length(hist_vals_a) >1
            smooth_data_g = 1.05*(smoothdata(interp1(hist_vals_a, 1:0.01:length(hist_vals_a))));
            fill((var_ind-1)*offset+[0 smooth_data_g 0], ...
                rescale((1:length([0 smooth_data_g 0]))/length([0 smooth_data_g 0]),...
                min(grp_a) - (bin_wd/10), max(grp_a) + (bin_wd/10)),[0.3010 0.7450 0.9330], 'FaceAlpha',0.3,...
                'LineStyle','none')
        end 
        plot((var_ind-1)*offset+[2,0], median(grp_a)*[1,1], "LineWidth", 2.5, "Color",  0.75*[0.3010 0.7450 0.9330] )
        
        if length(hist_vals_b) > 1
            smooth_data_b = 1.05*(smoothdata(interp1(hist_vals_b, 1:0.01:length(hist_vals_b))));
            fill((var_ind-1)*offset -[0 smooth_data_b 0],...
                rescale((1:length([0 smooth_data_b 0]))/length([0 smooth_data_b 0]),...
                min(grp_b) - (bin_wd/10), max(grp_b) + (bin_wd/10)), [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,...
                'LineStyle','none')
        end
        plot((var_ind-1)*offset+[-2,0], median(grp_b)*[1,1], "LineWidth", 2.5, "Color", 0.75*[0.8500 0.3250 0.0980] )
        set(gca, "XTick", [])
        if all(~isnan(y_lim))
            ylim(y_lim)
        end
        box off
        title(strrep(var, "_", " "))

        var_ind = var_ind +1;
    end

    nexttile
    hold on
    scatter(1,1,[],[0.3010 0.7450 0.9330], 'filled')
    scatter(1,0,[],[0.8500 0.3250 0.0980], 'filled')
    text(2,1, grp_names(1))
    text(2,0, grp_names(2))
    xlim([0,3])
    ylim([-15,3])
    hold off 
    axis off

    %legend(opts.grp_names, "Box","off", "Location","northeastoutside")
    if save_fig == 1
        saveas(f, sprintf("%s.%s", save_loc, file_type))
    end
end