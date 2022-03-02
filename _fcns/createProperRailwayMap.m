function [proper_railway_map,proper_measurement_selector,old_track_ids] = createProperRailwayMap(railway_map,measurement_data_selector,track_ids)
% [proper_railway_map,proper_measurement_selector,old_track_ids] = createProperRailwayMap(railway_map,measurement_data_selector,track_ids)
%

%% Initilization

if (nargin < 3) || isempty(track_ids)
    track_ids = railway_map.track_start_points.ID;
end % if
    
v_max = 160;

track_start_points_selector = ismember(railway_map.track_start_points.ID,track_ids);
track_maps_selector = ismember(railway_map.track_maps.ID,track_ids);

proper_railway_map.topology = railway_map.topology(track_start_points_selector,track_start_points_selector);
proper_railway_map.track_start_points = railway_map.track_start_points(track_start_points_selector,:);
proper_railway_map.track_maps = railway_map.track_maps(track_maps_selector,:);

%% Calculations

% Determine type of unkown tracks _________________________________________

for element_index = 1:size(proper_railway_map.track_maps,1)
    
    track_element_i = proper_railway_map.track_maps(element_index,:).track_element;    
    
    % Dismiss tracks that can't be uniquely identified
    if ~isnan(track_element_i) || element_index < 2 || element_index > (size(proper_railway_map.track_maps,1)-1)
        if isnan(track_element_i) && ( (element_index < 2) || (element_index > (size(proper_railway_map.track_maps,1)-1)) )
            warning('createProperRailwayMap: Dismissed a track-element! It can''t uniquely be identified!');
        end % of
        
        continue
    end % if
    
    % Determine track type by looking to the track's neighbors
    track_element_i_prev = proper_railway_map.track_maps(element_index-1,:).track_element;
    track_element_i_suc = proper_railway_map.track_maps(element_index+1,:).track_element;
    switch sprintf('%u%u',track_element_i_prev,track_element_i_suc)
        case '11'
            proper_railway_map.track_maps(element_index,:).track_element = 11;
            proper_railway_map.track_maps(element_index,:).r_start = 0;
            proper_railway_map.track_maps(element_index,:).r_end = 0;
            proper_railway_map.track_maps(element_index,:).speed_limit = calcSpeedLimit(v_max,0);
            
            p_0_prev = [proper_railway_map.track_start_points(element_index-1,:).x_0; ... 
                        proper_railway_map.track_start_points(element_index-1,:).y_0];
            phi_0_prev = proper_railway_map.track_start_points(element_index-1,:).phi_0;
            t_0_prev = [cosd(phi_0_prev); ... 
                        sind(phi_0_prev)];
            l_prev = proper_railway_map.track_maps(element_index-1,:).length;            
            p_end_prev = p_0_prev + l_prev*t_0_prev;          
                        
            x_0 = p_end_prev(1);
            y_0 = p_end_prev(2);
            x_end = proper_railway_map.track_start_points(element_index+1,:).x_0;
            y_end = proper_railway_map.track_start_points(element_index+1,:).y_0;
            phi_0 = atan2d(y_end-y_0,x_end-x_0);
            
            proper_railway_map.track_start_points(element_index,:).x_0 = x_0;
            proper_railway_map.track_start_points(element_index,:).y_0 = y_0;           
            proper_railway_map.track_start_points(element_index,:).phi_0 = phi_0;
            
        case '13'
            proper_railway_map.track_maps(element_index,:).track_element = 2;
            proper_railway_map.track_maps(element_index,:).r_start = 0;
            proper_railway_map.track_maps(element_index,:).r_end = proper_railway_map.track_maps(element_index+1,:).r_end;
            proper_railway_map.track_maps(element_index,:).speed_limit = calcSpeedLimit(v_max,proper_railway_map.track_maps(element_index,:).r_end);
        case '31'
            proper_railway_map.track_maps(element_index,:).track_element = 4;
            proper_railway_map.track_maps(element_index,:).r_start = proper_railway_map.track_maps(element_index-1,:).r_end;
            proper_railway_map.track_maps(element_index,:).r_end = 0;
            proper_railway_map.track_maps(element_index,:).speed_limit = calcSpeedLimit(v_max,proper_railway_map.track_maps(element_index,:).r_end);
        case '33'
            proper_railway_map.track_maps(element_index,:).track_element = 5;
            proper_railway_map.track_maps(element_index,:).r_start = proper_railway_map.track_maps(element_index-1,:).r_end;
            proper_railway_map.track_maps(element_index,:).r_end = proper_railway_map.track_maps(element_index+1,:).r_start;
            proper_railway_map.track_maps(element_index,:).speed_limit = calcSpeedLimit(v_max,proper_railway_map.track_maps(element_index,:).r_end);
    end % swich
    
end % for i

% Adjust measurement data selector for optimization _______________________

temp_selector = ismember(measurement_data_selector(:,1),track_ids);
proper_measurement_selector = measurement_data_selector(temp_selector,:);

% Discard unkown tracks at the start and the end __________________________

identified_tracks_tm_selector = ~isnan(proper_railway_map.track_maps.track_element');
identified_tracks_ids = unique(proper_railway_map.track_maps(identified_tracks_tm_selector,:).ID);
identified_tracks_sp_selector = ismember(proper_railway_map.track_start_points.ID,identified_tracks_ids);

proper_railway_map.topology = proper_railway_map.topology(identified_tracks_sp_selector,identified_tracks_sp_selector);
proper_railway_map.track_start_points = proper_railway_map.track_start_points(identified_tracks_sp_selector,:);
proper_measurement_selector = proper_measurement_selector(identified_tracks_sp_selector,:);
proper_railway_map.track_maps = proper_railway_map.track_maps(identified_tracks_tm_selector,:);

% Reset track IDs to be monotonically increasing __________________________

old_track_ids = identified_tracks_ids;
temp_map = proper_railway_map;
temp_z_selector = proper_measurement_selector;
for i = 1:length(old_track_ids)
    old_track_id = old_track_ids(i);
    new_track_id = i;
    
    sp_id_selector = ismember(proper_railway_map.track_start_points.ID,old_track_id);
    tm_id_selector = ismember(proper_railway_map.track_maps.ID,old_track_id);
    
    temp_map.track_start_points(sp_id_selector,:).ID = new_track_id;
    temp_map.track_maps(tm_id_selector,:).ID = new_track_id;
    
    temp_z_selector(sp_id_selector,1) = new_track_id;
end % for i
proper_railway_map = temp_map;
proper_measurement_selector = temp_z_selector;

end % function

