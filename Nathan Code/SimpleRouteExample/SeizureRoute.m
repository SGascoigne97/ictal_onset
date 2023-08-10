classdef SeizureRoute < handle
    properties
        touching_groups
        n_touching_groups
        n_states
        times
        wm_connections_used
        implanted
        touching
        wm_connections
        touching_without_implanted
        touching_using_unimplanted
        region_connections % 0 = no activation, 1 = initial region, 2 = touching, 
        % 3 = wm implanted-implanted, 4 = wm implanted-touching, 5 = wm touching-implanted, 6 = wm touching-touching
        % 7 = no connection

    end

    methods
        function obj = SeizureRoute(n_regions, n_states_total, implanted, touching, wm_connections)
            %SeizureRoute Construct an instance of this class
            %   initialises all the variables assuming the maximum number
            %   of possible groups is the number of implanted regionswm_connections = wm_connections;
            obj.touching_without_implanted = obj.touching;
            obj.touching_without_implanted(implanted,implanted) = false;
            obj.touching_using_unimplanted = (obj.distances(obj.touching_without_implanted) < 3) | obj.touching;
            obj.region_connections = zeros(n_regions, 1);
        end

        function newRegions(obj,r,t)
            %newRegions add new regions to the route
            % order agnostic
            % groups new regions by if they're touching
            % we will loop over every region adding it to all touching
            % groups that touch it
            % if multiple groups touch the region it will be added into
            % both groups

            obj.n_states = obj.n_states + 1;
            obj.times(obj.n_states) = t;
            bins = conncomp(graph(obj.touching_using_unimplanted(r,r)));
            added_touching = false(max(bins),1);
            for i_bin = 1:max(bins)
                % loop over all current groups
                for i_group = 1:obj.n_touching_groups
                    % if any of the regions in the group physically touch the
                    % new region, add it to that group
                    if any(obj.touching_using_unimplanted(obj.touching_groups(:,i_group,obj.n_states),r(bins == i_bin)))
                        obj.touching_groups(r(bins == i_bin),i_group,obj.n_states:end) = true;
                        added_touching(i_bin) = true;
                        obj.region_connections(r(bins == i_bin)) = 2;
                    end
                end
            end

            for i_bin = 1:max(bins)
                % if the region didn't touch any of the current groups
                if added_touching(i_bin)
                    continue
                end
                % create a new group with this region in
                new_regions = r(bins == i_bin);
                regions_in_groups = any(obj.touching_groups(:,:,obj.n_states),2);

                obj.touching_groups(new_regions, obj.n_touching_groups + 1, obj.n_states:end) = true;

                % mark any white matter connections to this region 
                regions_in_groups_plus_unimplanted_touching_regions = any(obj.touching_without_implanted(:,regions_in_groups),2) | regions_in_groups;
                unimplanted_touching_regions = find(any(obj.touching_without_implanted(new_regions,:),1));
                new_regions_plus_unimplanted_touching_regions = [new_regions', unimplanted_touching_regions];
                
                if any(obj.wm_connections(new_regions,regions_in_groups),'all')
                    % if there's WM from touching groups to new group
                    obj.wm_connections_used(new_regions,regions_in_groups,obj.n_states:end)...
                        = repmat(obj.wm_connections(new_regions,regions_in_groups),1,1,size(obj.wm_connections_used,3)-obj.n_states+1);
                    obj.region_connections(r(bins == i_bin)) = 3;
                elseif any(obj.wm_connections(new_regions,regions_in_groups_plus_unimplanted_touching_regions),"all")
                     % if there's WM from touching groups + adjacent unimplanted to new group
                    obj.wm_connections_used(new_regions,regions_in_groups_plus_unimplanted_touching_regions,obj.n_states:end)...
                        = repmat(obj.wm_connections(new_regions,regions_in_groups_plus_unimplanted_touching_regions),1,1,size(obj.wm_connections_used,3)-obj.n_states+1);
                    obj.region_connections(r(bins == i_bin)) = 4;
                elseif  any(obj.wm_connections(new_regions_plus_unimplanted_touching_regions,regions_in_groups),'all')
                    % if there's WM from touching groups to new group + adjacent unimplanted
                    obj.wm_connections_used(new_regions_plus_unimplanted_touching_regions,regions_in_groups,obj.n_states:end)...
                        = repmat(obj.wm_connections(new_regions_plus_unimplanted_touching_regions,regions_in_groups),1,1,size(obj.wm_connections_used,3)-obj.n_states+1);
                    obj.region_connections(r(bins == i_bin)) = 5;
                elseif any(obj.wm_connections(new_regions_plus_unimplanted_touching_regions,regions_in_groups_plus_unimplanted_touching_regions),"all")
                    % if there's WM + adjacent unimplanted from touching groups to new group + adjacent unimplanted
                    obj.wm_connections_used(new_regions_plus_unimplanted_touching_regions,regions_in_groups_plus_unimplanted_touching_regions,obj.n_states:end)...
                        = repmat(obj.wm_connections(new_regions_plus_unimplanted_touching_regions,regions_in_groups_plus_unimplanted_touching_regions),1,1,size(obj.wm_connections_used,3)-obj.n_states+1);
                    obj.region_connections(r(bins == i_bin)) = 6;
                else
                    obj.region_connections(r(bins == i_bin)) = 7;
                end


                obj.n_touching_groups = obj.n_touching_groups + 1;
                
                % this assigns a random touching group as the first if the
                % first state has multiple but it does still highlight the
                % precense of completely disconnected groups
               
                if obj.n_touching_groups == 1
                    obj.region_connections(r(bins == i_bin)) = 1;
                end
            end
        end
        function simpleBrainPlot(obj, atlas, atlasToUse, t, crop) 
            roi_names = atlas.name{atlasToUse};
            %roi_names = strrep(roi_names, 'r.', 'ctx-rh-');
            %roi_names = strrep(roi_names, 'l.', 'ctx-lh-');
            values = zeros(size(roi_names));
            for i_group = 1:obj.n_touching_groups
                values(obj.touching_groups(:,i_group,t)) = i_group;
            end
            brewermap('Set1');
            colors = brewermap(max(values));
            plotBrain_png(roi_names(values > 0),values(values > 0),colors);
        end
        function visualise(obj, atlasToUse, atlas, lf, lv, rf, rv)
        
            xyz = atlas.xyz{atlasToUse,1};
            fig = figure('units','pixels','Position',[2151 141 941 700]);
            sld = uicontrol('Style', 'slider',...
                'Min',1,'Max',obj.n_states,'Value',1,...
                'Units', 'Normalized',...
                'Position', gca().Position + [-0.02 -0.87 0.04 0],...
                'SliderStep', [1/(obj.n_states-1+2) 0.01],...
                'Callback', @value_changed);
            setup = false;
            plot_touching_groups(1);
            obj.simpleBrainPlot(atlas, atlasToUse, 1, 0);

            function plot_touching_groups(t)
                subplot(1,2,1)
                if setup
                    pos = fig.CurrentAxes.CameraPosition;
                    target = fig.CurrentAxes.CameraTarget;
                end
                scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'black');
                hold on
                scatter3(xyz(obj.implanted,1),xyz(obj.implanted,2),xyz(obj.implanted,3),'black','filled');
                
                brewermap('Set1');
                colors = brewermap(obj.n_touching_groups);
                for i_group = 1:obj.n_touching_groups
                    group = obj.touching_groups(:,i_group,t);
                    scatter3(xyz(group,1), xyz(group,2), xyz(group,3),30,colors(i_group, :), 'filled');
                    %hold on
                    for i = find(group)'
                        for j = find(group)'
                            if obj.touching(i,j)
                                x = [xyz(i,1), xyz(j,1)];
                                y = [xyz(i,2), xyz(j,2)];
                                z = [xyz(i,3), xyz(j,3)];
                                h = line(x,y,z,'color',colors(i_group, :),'LineWidth',1);
                            end
                        end
                    end
                end
                [x,y] = find(obj.wm_connections_used(:,:,t));
                for i = 1:length(x)
                    x_coord = [xyz(x(i),1), xyz(y(i),1)];
                    y_coord = [xyz(x(i),2), xyz(y(i),2)];
                    z_coord = [xyz(x(i),3), xyz(y(i),3)];
                    h = line(x_coord,y_coord,z_coord,'color','black','LineWidth',1);
                end
                Hl = patch('vertices',lv,'faces',lf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
                Hr = patch('vertices',rv,'faces',rf(:,[1 3 2]),'facecolor',[.5 .5 .5],'facealpha',0.1,'edgecolor','none');
                camlight('headlight','infinite')
                axis equal
                hold off
                if setup
                    fig.CurrentAxes.CameraPosition = pos;
                    fig.CurrentAxes.CameraTarget = target;
                end
                setup = true;
                title('t = ' + string(obj.times(t) - obj.times(1)))
            end
            function value_changed(hObject,callbackdata)
                newval = sld.Value;                         %get value from the slider
                newval = round(newval);                         %round off this value
                if newval ~= sld.Value
                    set(sld, 'Value', newval);                  %set slider position to rounded off values
                end
                plot_touching_groups(newval);
                subplot(1,2,2)
                obj.simpleBrainPlot(atlas, atlasToUse, newval, 0);
            end

        end
    end
    methods(Access=protected)
        function distances = distances(obj, touching)
            distances = inf(length(touching));
            distances(touching == 1) = 1;
            for v = 1:length(touching)
                distances(v,v) = 0;
            end
            for k = 1:length(touching)
                for i = 1:length(touching)
                    for j = 1:length(touching)
                        if distances(i,j) > distances(i,k) + distances(k,j) 
                            distances(i,j) = distances(i,k) + distances(k,j);
                        end
                    end
                end
            end
        end

    end
end