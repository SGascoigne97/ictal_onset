pat = 3;
sz = 2;
ch_2_roi = final_output.chan_2_roi_matrix_120{pat,1};
roi = repmat(final_output.imprint_roi_120{pat,1}(:,sz),1,size(final_output.chan_2_roi_matrix_120{pat,1},2));
chan = repmat(final_output.imprint_chan{pat,1}(:,sz)', size(final_output.chan_2_roi_matrix_120{pat,1},1),1);

plot_mat = 3*ch_2_roi+roi+chan; 
plot_mat(plot_mat == 2) = 1;
plot_mat(plot_mat == 4) = 3;
plot_mat(plot_mat == 5) = 4;

figure(3)
subplot(5,5, [1:4,6:9,11:14,16:19])
imagesc(plot_mat)
clim([0 4])
set(gca, "XTick", [])
ylabel("Regions (Laus 120)")
subplot(5,5,5:5:20)
imagesc(final_output.imprint_roi_120{pat,1}(:,sz))
clim([0,4])
set(gca, "YTick", [], "XTick", [])
subplot(5,5,21:24)
imagesc(final_output.imprint_chan{pat,1}(:,sz)')
set(gca, "YTick", [])
clim([0,4])
xlabel("Channels")
sgtitle(sprintf("%s Seizure %d",final_output.Patient_id{pat,1}, sz ))

%%
figure()
subplot(121)
imagesc(final_output.clo_roi_120{pat,1})
subplot(122)
imagesc(final_output.imprint_roi_120{pat,1}(:,sz))
%%
sum(final_output.clo_roi_120{pat,1})
sum(final_output.imprint_roi_120{pat,1}(:,sz))
