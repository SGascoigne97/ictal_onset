function seizureRoute = showRoute(wm_names, region_imprint, touching_wm_size, wm_touching, region_indexes, roiIndex, atlas, lf, lv, rf, rv)
    %seizureRoute = SeizureRoute(size(wm_names,1),size(region_imprint_no_duplicates,1),size(region_imprint_no_duplicates,2),touching_wm_size,wm_touching);
    seizureRoute = SeizureRoute(size(wm_names,1), size(region_imprint,1), region_indexes, touching_wm_size, wm_touching);
    for i = 1:size(region_imprint,2)
        if i > 1
            regions = (region_imprint(:,i) - any(region_imprint(:,1:i-1),2)) > 0;
        else 
            regions = region_imprint(:,i);
        end
        if any(regions)
            seizureRoute.newRegions(region_indexes(regions),i/8);
            
        end
    end
    seizureRoute.visualise(roiIndex, atlas, lf, lv, rf, rv);
end
