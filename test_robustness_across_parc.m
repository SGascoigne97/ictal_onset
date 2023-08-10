
onset_across_titles = ["per_sz", "across_sz"];
det_meths = ["imprint", "EI"];
for onset_across = 1%[0,1]
    for summ_meas = ["med", "max"]
        for det_meth = det_meths
            fig = figure("Position",[10,10,700,900]);
            if onset_across == 0
                roi_120 = comp_meth_with_clo.(sprintf(onset_across_titles(onset_across+1))).roi_120.(sprintf(det_meth)).(sprintf(summ_meas));
                roi_250 = comp_meth_with_clo.(sprintf(onset_across_titles(onset_across+1))).roi_250.(sprintf(det_meth)).(sprintf(summ_meas));
            else
                roi_120 = comp_meth_with_clo.(sprintf(onset_across_titles(onset_across+1))).roi_120.(sprintf(det_meth));
                roi_250 = comp_meth_with_clo.(sprintf(onset_across_titles(onset_across+1))).roi_250.(sprintf(det_meth));
            end
            roi_comp_tab = join(roi_120,roi_250, 'Keys', "Patient_id");
            roi_comp_tab = roi_comp_tab(~isnan(roi_comp_tab.Jacc_z_roi_120),:);
            roi_comp_tab = roi_comp_tab(~isnan(roi_comp_tab.Jacc_z_roi_250),:);
            roi_comp_tab = roi_comp_tab(~isnan(roi_comp_tab.Perc_z_roi_250),:);
            
            for comp_measure = comp_measures
                nexttile
                val_120 = roi_comp_tab.(sprintf("%s_roi_120", comp_measure));
                val_250 = roi_comp_tab.(sprintf("%s_roi_250", comp_measure));
%                 [rho,p] = corr(val_120,val_250, 'type','Spearman');
                
                % Compute the distance between 250 and 120 to capture the sum of squared errors 
                err = val_250 - val_120;
                sse = sqrt(sum(err.^2)/length(err));
                scaled_sse = sse/mean(val_120);

                scatter(val_120, val_250, 'filled')
                xlim([min([val_120;val_250]), max([val_120;val_250])])
                ylim([min([val_120;val_250]), max([val_120;val_250])])
                axis square
               
                hold on
                plot([min(val_120), max(val_120)],[min(val_120), max(val_120)])
                hold off
                xlabel("Lausanne 120")
                xlabel("Lausanne 250")
                if contains(comp_measure, "z")
                    patch([-2 2 2 -2],[-2 -2 2 2],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
                    text(0,1.5,"chance", "HorizontalAlignment","center")
                end
%                 title(sprintf("%s (rho = %.2f, p = %.5f)", comp_measure, rho, p))
                title(sprintf("%s (sum sq error %.3f)", comp_measure, scaled_sse))
            end
            if onset_across == 0
                sgtitle(sprintf("%s subject %s across parcellations", det_meth, summ_meas))
                saveas(fig,sprintf("figures/rapid_prototype/robust/%s_subject_%s.png", det_meth, summ_meas))
            else
                sgtitle(sprintf("%s comparing across seizure onset across parcellations", det_meth))
                saveas(fig,sprintf("figures/rapid_prototype/robust/%s_onset_across.png", det_meth))
            end
        end
    end
end