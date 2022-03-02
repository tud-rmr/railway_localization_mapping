% 
% Prepare data for map creation via optimization
%   
%   Author: Hanno Winter
%   Date: 07-Apr-2021; Last revision: 09-May-2021

%% Init

if ~exist('init_localization_executed','var')
    error('prepareOptimInData.m: Initaliziation script not executed!')
end % if
if ~exist('simout_optimization_data_selector','var')
    error('prepareOptimInData.m: SimOut data not available!')
end % if

fprintf('Preparing input data for map creation via optimization')

%% Load data and prepare it for optimization

% Load measurement data ___________________________________________________

switch optim_input_selector
    case {'GPS'}       
        optim_input_data = [sim_gnss_filter_input.UtmEast_m';sim_gnss_filter_input.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
        optim_input_data_P = cat(3,sim_gnss_filter_input.PositionCov{:});
        optimization_data_selector = optimization_gps_data_selector;
    case {'IMM'}        
        optim_input_data = [sim_imm.UtmEast_m';sim_imm.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
        optim_input_data_P = cat(3,sim_imm.PositionCov{:});
        optimization_data_selector = optimization_imm_data_selector;
	case {'TGC'}        
        optim_input_data = [sim_tgc.UtmEast_m';sim_tgc.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
        optim_input_data_P = cat(3,sim_tgc.PositionCov{:});
        optimization_data_selector = optimization_imm_data_selector;
    otherwise
        error('map_optimization: Input data string not supported!');
end % switch

% Create proper map from simulation data __________________________________

[proper_sim_map,optim_data_selector,old_track_ids] = ... 
    createProperRailwayMap(sim_map,optimization_data_selector,optim_track_ids);

% Reduce measurment points ________________________________________________

selection_mask = ismember((1:size(optim_data_selector,2)-1),(1:z_data_dividor:size(optim_data_selector,2)-1));
for i = 1:size(optim_data_selector,1)
    optim_data_selector(i,2:end) = optim_data_selector(i,2:end) & selection_mask;
end % for i

%% Determine optimization frames

% Select only straights before and after circular arcs ____________________

% Select straights and arcs
% straights_selector = (proper_sim_map.track_maps.track_element==1);
straights_selector = (proper_sim_map.track_maps.track_element==1) | (proper_sim_map.track_maps.track_element==11);
arc_selector = (proper_sim_map.track_maps.track_element==3);
element_ids = proper_sim_map.track_maps.ID((straights_selector | arc_selector),:);
old_element_ids = old_track_ids((straights_selector | arc_selector));

% Find standstill for each element (on original railway-map)
v_zero_selector = (sim_imm.VelocityVehicle_ms == 0);
element_dir_temp = zeros(length(proper_sim_map.track_maps.ID),1);
for i = 1:length(proper_sim_map.track_maps.ID)
    %old_element_id_i = proper_sim_map.track_maps(i,:).ID;
    old_element_id_i = old_track_ids(i);
    z_i_selector = (optimization_imm_data_selector(:,1) == old_element_id_i);
    z_i = optimization_imm_data_selector(z_i_selector,2:end);
    
    if sum(z_i(v_zero_selector)) > 0 
        element_dir_temp(i) = 1;
    end % if
end % for

% Find elements with standstill on the new proper railway-map, which is used for optimization
% old_stillstand_ids = proper_sim_map.track_maps.ID(logical(element_dir_temp));
old_stillstand_ids = old_track_ids(logical(element_dir_temp));
element_dir = zeros(length(element_ids),1);
for i = 1:sum(element_dir_temp)
    old_id = old_stillstand_ids(i);
    new_id_idx = find(old_element_ids==old_id,1,'first');   
    
    while (old_id > 0) && isempty(new_id_idx)
        old_id = old_id -1;
        new_id_idx = find(old_element_ids==old_id,1,'first'); 
    end % while    

    if ~isempty(new_id_idx)
        element_dir(new_id_idx) = 1;        
    end % if
end % for i

% Set length to zero for elements with standstill
if any(element_dir)
    standstill_ids = element_ids(logical(element_dir));
    
    standstill_selector = ismember(proper_sim_map.track_maps.ID,standstill_ids);
    proper_sim_map.track_maps.length(standstill_selector) = 0;
end % if

% Save indices where switch of traveling or standstill occurs
drive_end_indices = sort( ... 
                          [strfind(num2str(element_dir)','1')'; ... 
                           length(element_ids)], ... 
                          'asc' ... 
                        );
                    
% Create list of track IDs usable for optimization (and thereby respect the driving direction)
optim_indices = cell(length(drive_end_indices),1);
optim_frame_size = optim_indices;
optim_start_indices = optim_indices;
optim_end_indices = optim_indices;
package_start_indices = optim_indices;
package_end_indices = optim_indices;
for drive_i = 1:length(drive_end_indices)
    
    if drive_i == 1
        drive_indices = 1:drive_end_indices(drive_i);
    else
        drive_indices = (drive_end_indices(drive_i-1)+1):drive_end_indices(drive_i);
    end % if    
    
    element_ids_drive_i = element_ids(drive_indices);        
    track_map_h = [];
    for i = 1:length(element_ids_drive_i)

        id_i = element_ids_drive_i(i);
        index_i = find(proper_sim_map.track_maps.ID == id_i);
        track_map_i = proper_sim_map.track_maps(index_i,:);

        if (track_map_i.track_element == 1) ... % Current: straight; Previous: arc
           && ... 
           (isempty(track_map_h) || (track_map_h.track_element == 3))

            optim_indices{drive_i}(end+1,1) = index_i;

        elseif (track_map_i.track_element == 1) ... % Current: straight; Previous: straight
                && ... 
                (track_map_h.track_element == 1) % ... 
%                 && ... 
%                 (track_map_i.length > proper_sim_map.track_maps(optim_indices{drive_i}(end,1),:).length)

            optim_indices{drive_i}(end+1,1) = index_i;
            
        elseif (track_map_i.track_element == 3) ... % Current: Arc; Previous: straight
           && ... 
           (isempty(track_map_h) || (track_map_h.track_element == 1))
       
            optim_indices{drive_i}(end+1,1) = index_i;
       
        elseif (track_map_i.track_element == 3) ... % Current: Arc; Previous: Arc
                && ... 
                (track_map_h.track_element == 3) % ... 
%                 && ... 
%                 (track_map_i.length > proper_sim_map.track_maps(optim_indices{drive_i}(end,1),:).length)
            
            optim_indices{drive_i}(end+1,1) = index_i;

        end % if
        
        if track_map_i.track_element ~= 11
            track_map_h = track_map_i;
        end % if

    end % for i

    % Select only straights were measurements are available ___________________

    poor_z_selector = sum(optim_data_selector(optim_indices{drive_i},2:end),2) < min_z;
    optim_indices{drive_i} = optim_indices{drive_i}(~poor_z_selector);

    % Create optimization frames ______________________________________________

    optim_frame_size{drive_i} = min(frame_size,length(optim_indices{drive_i}));
    optim_start_indices{drive_i} = optim_indices{drive_i}(1:end-(optim_frame_size{drive_i}-1));
    optim_end_indices{drive_i} = optim_indices{drive_i}(optim_frame_size{drive_i}:end);
    
    % Fade in optimization frame
    if 0
        optim_start_indices{drive_i} = [optim_indices{drive_i}(1:optim_frame_size{drive_i}-2);optim_start_indices{drive_i}];
        optim_end_indices{drive_i} = [optim_indices{drive_i}(3:optim_frame_size{drive_i});optim_end_indices{drive_i}];
    end % if
    
    % Fade out optimization frame
    if 0
        optim_start_indices{drive_i} = [optim_start_indices{drive_i};optim_indices{drive_i}(end-optim_frame_size{drive_i}+2:end-1)];
        optim_end_indices{drive_i} = [optim_end_indices{drive_i};repmat(optim_end_indices{drive_i}(end),optim_frame_size{drive_i}-2,1)];
    end % if
    
    % Create output map frames ____________________________________________

    package_start_indices{drive_i} = optim_start_indices{drive_i};
    package_end_indices{drive_i} = optim_end_indices{drive_i};
    
    % Add one more element to package
    if 1
        package_start_indices{drive_i}(2:end) = package_start_indices{drive_i}(2:end)+1;
        package_end_indices{drive_i} = package_end_indices{drive_i};    
    end % if
        
end % for

% Convert per drive list to single list
optim_indices = cell2mat(optim_indices);
optim_start_indices = cell2mat(optim_start_indices);
optim_end_indices = cell2mat(optim_end_indices);
package_start_indices = cell2mat(package_start_indices);
package_end_indices = cell2mat(package_end_indices);

% Trim input map __________________________________________________________

optim_map_indices = (optim_indices(1):optim_indices(end));
optim_map.topology = proper_sim_map.topology(optim_map_indices,optim_map_indices);
optim_map.track_start_points = proper_sim_map.track_start_points(optim_map_indices,:);
optim_map.track_maps = proper_sim_map.track_maps(optim_map_indices,:);
optim_data_selector = optim_data_selector(optim_map_indices,:);

% Adjust frame indices to trimmed map _____________________________________

idx_offset = optim_start_indices(1)-1;
optim_start_indices = optim_start_indices-idx_offset;
optim_end_indices = optim_end_indices-idx_offset;
package_start_indices = package_start_indices-idx_offset;
package_end_indices = package_end_indices-idx_offset;

% Final optimization run over full map
if 0
    optim_start_indices = [optim_start_indices;1];
    optim_end_indices = [optim_end_indices;optim_end_indices(end)];
    package_start_indices = [package_start_indices;optim_start_indices(end)];
    package_end_indices = [package_end_indices;optim_end_indices(end)];
end % if

%% Finish

% nothing to do
