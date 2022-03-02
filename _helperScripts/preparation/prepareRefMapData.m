% 
% Prepare reference map data
%   
%   Author: Hanno Winter
%   Date: 11-Apr-2021; Last revision: 18-Jul-2021

%% Settings

v_cdf_min = v_min_tgc_speed-v_min_tgc_hysteresis; % min speed to use for CDF calculations

%% Init

if ~exist('map_data_prepared','var') || map_data_prepared == false
    map_data_prepared = false;
    fprintf('Preparing map data');
elseif map_data_prepared
    fprintf('Already prepared map data\n');
    return  
end % if

prepareSimOutData

%% Calculations

switch input_data_selector
    
    case {'C'}        
        
        importCMapData
        
        % Calculate track-map selector to limit it to the availabe
        % measurements ____________________________________________________
        
        ref_track_map_selector = limitRefMap(ref_track_map,sim_imm);
        ref_track_map_indices = find(ref_track_map_selector);
        
        % Relate IMM data to reference track-map __________________________
        
        imm_ref_track_map = relateRefMapToImmData(ref_track_map,ref_track_map_selector,sim_imm);
 
        % Notify that reference track-map is available ____________________
        
        plot_ref_map = true;
        
    case {'BS'}
        
        importBsMapData
        
        % Manually set list of relevant track IDs _________________________
        
        relevant_track_ids = [7 23 24 21 11 17 16 15 10 25];
        relevant_track_ids = fliplr(relevant_track_ids)';
        
        % Create reference track-map ______________________________________
        
        % Create reference trac-map table
        ref_track_map = table();
        for i = 1:length(relevant_track_ids)
            
            track_id_i = relevant_track_ids(i);
            
            % Shrink map
            track_map_selector = ismember(map_xy_track_map.track_num,track_id_i);

            % Write to table
            ref_track_map_data = table();
            ref_track_map_data.Latitude_deg = [map_xy_track_map(track_map_selector,:).lat];
            ref_track_map_data.Longitude_deg = [map_xy_track_map(track_map_selector,:).lon];
            ref_track_map_data.UtmEast_m = [map_xy_track_map(track_map_selector,:).utm_east];
            ref_track_map_data.UtmNorth_m = [map_xy_track_map(track_map_selector,:).utm_north];
            
            ref_track_map = [ref_track_map;ref_track_map_data];
        end % if     
                        
        % Calc reference track-map distance 
        delta_d = [[0;0],diff([ref_track_map.UtmEast_m,ref_track_map.UtmNorth_m]',1,2)];
        ref_track_map.TrackDistance_m = cumsum(sqrt(sum(delta_d.^2,1)))';
        
        % Limit map to measurement data 
        ref_track_map_selector = true(size(ref_track_map,1),1);
        
        % Relate IMM data to reference track-map __________________________ 
        
        % imm_ref_track_map = relateRefMapToImmData(ref_track_map,ref_track_map_selector,sim_imm);
        imm_ref_track_map = [];
        
        % Notify that reference track-map is available ____________________
        plot_ref_map = true;
        
        
    otherwise % no reference track-map availalbe
        
        ref_track_map = [];
        ref_track_map_selector = [];
        imm_ref_track_map = [];
        plot_ref_map = false;
        
end % switch

%% Deviation from maps

calcDeviationToRefMap

%% Finish

map_data_prepared = true;

%% Save

prompt_str = 'Do you want to save the processed map data (j/n)?';
user_str = input(prompt_str,'s');

switch user_str
    case {'j','y','yes','ja'}
        map_data_saved = false;
        save_map_data = true;
    case {'n','no','nein'}
        save_map_data = false;
    otherwise
        save_map_data = false;
end % switch

if save_map_data
    saveMapData
end % if

%% Helper functions

function ref_map_selector = limitRefMap(ref_track_map,pos_data)

[~,min_lat_idx] = min(pos_data.Latitude_deg);
[~,max_lat_idx] = max(pos_data.Latitude_deg);
[~,min_lon_idx] = min(pos_data.Longitude_deg);
[~,max_lon_idx] = max(pos_data.Longitude_deg);

p_start_idx = min([min_lat_idx,max_lat_idx,min_lon_idx,max_lon_idx]);
p_end_idx = max([min_lat_idx,max_lat_idx,min_lon_idx,max_lon_idx]);

p_start = [pos_data.Latitude_deg(p_start_idx);pos_data.Longitude_deg(p_start_idx)];
p_end = [pos_data.Latitude_deg(p_end_idx);pos_data.Longitude_deg(p_end_idx)];
p_ref = [ref_track_map.Latitude_deg';ref_track_map.Longitude_deg'];

p_start_diff = p_ref - p_start;
p_start_diff = sqrt(p_start_diff(1,:).^2+p_start_diff(2,:).^2);
[~,start_idx] = min(p_start_diff);

p_end_diff = p_ref - p_end;
p_end_diff = sqrt(p_end_diff(1,:).^2+p_end_diff(2,:).^2);
[~,end_idx] = min(p_end_diff);

ref_map_selector = ismember(1:size(p_ref,2)',min(start_idx,end_idx):max(start_idx,end_idx));

end % function

function imm_ref_track_map = relateRefMapToImmData(ref_track_map,ref_track_map_selector,sim_imm)
    
fprintf('Relate IMM data to reference track-map:\n');

% Settings
d_x = 500; % search radius for corresponding reference points in m
d_y = 500; % search radius for corresponding reference points in m

% Find nearest points on reference-map
p_ref = [ref_track_map.UtmEast_m';ref_track_map.UtmNorth_m'];
p_test = [sim_imm.UtmEast_m';sim_imm.UtmNorth_m'];
imm_track_length = nan(size(sim_imm.Time,1),1);     
for i = 1:size(p_test,2)
    p_test_i = p_test(:,i);

    % Status
    waitbarStatus(i-1,length(imm_track_length),5)

    % Pre-select reference data
    x_ref_data_selector = (p_ref(1,:) > p_test_i(1)-d_x) & (p_ref(1,:) < p_test_i(1)+d_x);
    y_ref_data_selector = (p_ref(2,:) > p_test_i(2)-d_y) & (p_ref(2,:) < p_test_i(2)+d_y);
    ref_data_selector = x_ref_data_selector & y_ref_data_selector;
    ref_data_indices = find(ref_data_selector);

    % Find nearest points on ref-track
    deltas = sqrt(sum((p_ref(:,ref_data_selector)-p_test_i).^2,1));
    [~,min_delta_idx] = min(deltas);            

    projection_on_ref_track = nan(sum(ref_data_selector),1);
    %for j = 1:sum(ref_data_selector)
    for j = [max(min_delta_idx-1,1),min_delta_idx,min(min_delta_idx+1,size(ref_data_indices,2))]
        ref_data_index_j = ref_data_indices(j);

        % p_ref_j = p_ref(:,j);                
        p_ref_j = p_ref(:,ref_data_index_j);
        if ref_data_index_j < size(p_ref,2)
            t_ref_j = p_ref(:,ref_data_index_j+1) - p_ref(:,ref_data_index_j);
        else
            t_ref_j = p_ref(:,ref_data_index_j) - p_ref(:,ref_data_index_j-1);
        end % if       

        projection_on_ref_track(j) = dot(p_test_i-p_ref_j,t_ref_j)/sum(t_ref_j.^2);
    end % for j           

    [~,lambda_min_idx] = min(abs(projection_on_ref_track));
    if isempty(lambda_min_idx)
        continue
    end
    lambda_min = projection_on_ref_track(lambda_min_idx);

    % Calculate track-length of mapped IMM point            
    ref_map_idx = ref_data_indices(lambda_min_idx);
    imm_track_length(i) = ref_track_map.TrackDistance_m(ref_map_idx) + lambda_min;

end % for i
waitbarStatus(i,length(imm_track_length),5)

% Store only unique data points in map 
[~,ia,~] = unique(ref_track_map.TrackDistance_m);
ref_track_map = ref_track_map(ia,:);
ref_track_map_selector = ref_track_map_selector(ia);

% Create map with reference to IMM data points
imm_ref_track_map = ... 
    interpTrackMap(imm_track_length,ref_track_map(ref_track_map_selector,:).TrackDistance_m,ref_track_map(ref_track_map_selector,:));

end % function
