% Create figure displaying what each comparison measure means
load('roi_info/lhSurf.mat') 
load('roi_info/rhSurf.mat')
%%
patients = string(final_comp.Patient_id);

for pat = 1:length(patients)
    patient = patients(pat);
    pat_onset = final_output(final_output.Patient_id == patient,:);
    
    % xyz coordinates for regions in 125 atlas
    
    xyz = atlas.xyz{3,1};
    roi_names = atlas.name{3,1};
    
    patient_rois = pat_onset.ROI_ids{:};
    incl_roi = ismember(roi_names, patient_rois);
    patient_xyz = xyz(incl_roi,:); % Careful - the order of regions is not preserved here 
    
    labelled_onset = pat_onset.labelled_onset_roi{1,1};
    labelled_onset(isnan(labelled_onset)) = 0;
    onset_names = patient_rois(find(labelled_onset));
    resected = pat_onset.resected_roi{1,1};
    resected(isnan(resected)) = 0;
    resected_names = patient_rois(find(resected));
    onset_and_resec = labelled_onset + resected == 2;
    o_and_r_names = patient_rois(find(onset_and_resec));
    
    figure(pat)
    set(gcf,'Position', [10 10 900 600])
    
    % Map out resection
    subplot(331)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, resected_names)),1),...
        xyz(find(ismember(roi_names, resected_names)),2),...
        xyz(find(ismember(roi_names, resected_names)),3), 'k','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    view([90 0])
    axis equal
    axis off
    title("Resected regions")
    subplot(332)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, resected_names)),1),...
        xyz(find(ismember(roi_names, resected_names)),2),...
        xyz(find(ismember(roi_names, resected_names)),3), 'k','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    view([0 90])
    axis equal
    axis off
    subplot(333)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, resected_names)),1),...
        xyz(find(ismember(roi_names, resected_names)),2),...
        xyz(find(ismember(roi_names, resected_names)),3), 'k','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    axis equal
    axis off
    view([270 0])
   
    % Map out CLO
    subplot(334)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, onset_names)),1),...
        xyz(find(ismember(roi_names, onset_names)),2),...
        xyz(find(ismember(roi_names, onset_names)),3), 'r','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    view([90 0])
    axis equal
    axis off
    title("Clinically labelled onset regions")
    subplot(335)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, onset_names)),1),...
        xyz(find(ismember(roi_names, onset_names)),2),...
        xyz(find(ismember(roi_names, onset_names)),3), 'r','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    view([0 90])
    axis equal
    axis off
    subplot(336)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, onset_names)),1),...
        xyz(find(ismember(roi_names, onset_names)),2),...
        xyz(find(ismember(roi_names, onset_names)),3), 'r','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    axis equal
    axis off
    view([270 0])
    
    % Map out onset and resection
    subplot(337)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, resected_names)),1),...
        xyz(find(ismember(roi_names, resected_names)),2),...
        xyz(find(ismember(roi_names, resected_names)),3), 'k','filled')
    scatter3(xyz(find(ismember(roi_names, onset_names)),1),...
        xyz(find(ismember(roi_names, onset_names)),2),...
        xyz(find(ismember(roi_names, onset_names)),3), 'r','filled')
    scatter3(xyz(find(ismember(roi_names, o_and_r_names)),1),...
        xyz(find(ismember(roi_names, o_and_r_names)),2),...
        xyz(find(ismember(roi_names, o_and_r_names)),3), 'b','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    view([90 0])
    axis equal
    axis off
    title("Onset and resection")
    subplot(338)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, resected_names)),1),...
        xyz(find(ismember(roi_names, resected_names)),2),...
        xyz(find(ismember(roi_names, resected_names)),3), 'k','filled')
    scatter3(xyz(find(ismember(roi_names, onset_names)),1),...
        xyz(find(ismember(roi_names, onset_names)),2),...
        xyz(find(ismember(roi_names, onset_names)),3), 'r','filled')
    scatter3(xyz(find(ismember(roi_names, o_and_r_names)),1),...
        xyz(find(ismember(roi_names, o_and_r_names)),2),...
        xyz(find(ismember(roi_names, o_and_r_names)),3), 'b','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    view([0 90])
    axis equal
    axis off
    subplot(339)
    scatter3(patient_xyz(:,1),patient_xyz(:,2),patient_xyz(:,3), 'k')
    hold on
    scatter3(xyz(find(ismember(roi_names, resected_names)),1),...
        xyz(find(ismember(roi_names, resected_names)),2),...
        xyz(find(ismember(roi_names, resected_names)),3), 'k','filled')
    scatter3(xyz(find(ismember(roi_names, onset_names)),1),...
        xyz(find(ismember(roi_names, onset_names)),2),...
        xyz(find(ismember(roi_names, onset_names)),3), 'r','filled')
    scatter3(xyz(find(ismember(roi_names, o_and_r_names)),1),...
        xyz(find(ismember(roi_names, o_and_r_names)),2),...
        xyz(find(ismember(roi_names, o_and_r_names)),3), 'b','filled')
    patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none')
    axis equal
    axis off
    view([270 0])
    
    sgtitle(sprintf("Patient %s", patient))
    vis_save_plot(gcf, '/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/figures/onset_resec/', ...
        sprintf('%s',patient), {'png'})
end


