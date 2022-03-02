function [topology,track_maps,track_start_points] = x2contMatRailwayMap(x,x_descriptor,start_point)
% [topology,track_maps,track_start_points] = x2contMatRailwayMap(x,x_descriptor,start_point)
% 

%% Settings 

v_max = 160; % Max speed on straights in km/h

%% Checks

if (size(x,1) ~= 4) 
    error('X2ContRailwayMap: Wrong parameter vector ''x''!');
end 

if sum(x_descriptor == 0) > 2
    error('X2ContRailwayMap: Too many railway-map start-points declared!');
end % if

if (nargin < 3)
    start_point = [];
end % if

if ~any(x_descriptor == 0) && (x_descriptor(1) ~= 1) && isempty(start_point)
    error('X2ContRailwayMap: No railway-map start-point declared!');
end % if

%% Init

straight_idx = 1;
transition_straight_idx = 11;
arc_idx = 3;

tm_id_idx = 1;
tm_element_idx = 2;
tm_r0_idx = 3;
tm_rend_idx = 4;
tm_length_idx = 5;
tm_vmax_idx = 6;

sp_id_idx = 1;
sp_x0_idx = 2;
sp_y0_idx = 3;
sp_phi0_idx = 4;

empty_track_start_point = @(id) [id nan nan nan nan];
empty_track_map = @(id) [id nan nan nan nan nan nan];

% Respect transition straights / fit x to length of x_descriptor
x = restoreX(x,x_descriptor);

%% Calculations I/II

% Adapt, if start or end-track are circular arcs (needed for optimizer to run smoothly)
if x_descriptor(1) == 0
    x(1,2) = 0;
end % if
if x_descriptor(end) == 99
    x(4,end) = 0;
end % if

% Extract track element order from railway-map vector _____________________

track_elements = [];
i = 1;
while i <= size(x,2)
    
    % Straight
    if (x_descriptor(1,i) == straight_idx) 
        
        track_elements = [track_elements;1];
        i = i + 1;
        
    elseif (x_descriptor(1,i) == transition_straight_idx)
        
        track_elements = [track_elements;11];
        i = i + 1;
    
	% Circular arc
    elseif (x_descriptor(1,i) == arc_idx) || (x_descriptor(1,i) == 13) || (x_descriptor(1,i) == 31) || (x_descriptor(1,i) == 131)
        
        % Force direct transition from straight to arc
        if x_descriptor(1,i) == 13 || (x_descriptor(1,i) == 131)
            x(1,i) = 0;
        end % if
        
        % Force direct transition from arc to straight
        if x_descriptor(1,i) == 31 || (x_descriptor(1,i) == 131)
            x(4,i) = 0;
        end % if
        
        % Set track-element order
        if ~( x(1,i) == 0 ) 
            track_elements = [track_elements;2;3];
        else
            track_elements = [track_elements;3];
        end % if
        
        % Look for succeeding arcs and fill with clothoids
        i = i+1;
        while (i <= size(x_descriptor,2)) && ( (x_descriptor(1,i) == arc_idx) || (x_descriptor(1,i) == 13) || (x_descriptor(1,i) == 31) )
            
            % Force direct transition from straight to arc
            if x_descriptor(1,i) == 13 || (x_descriptor(1,i) == 131)
                x(1,i) = 0;
            end % if

            % Force direct transition from arc to straight
            if x_descriptor(1,i) == 31 || (x_descriptor(1,i) == 131)
                x(4,i) = 0;
            end % if
            
            % Set track-element order
            track_elements = [track_elements;5;3];
            
            i = i+1;            
        end % while j
        
        % Finish series of arcs
        %if ~isnan(x(4,i-1))
        if ~(x(4,i-1) == 0)
            track_elements = [track_elements;4];
        end % if    
    
	% Start-point
    elseif (x_descriptor(1,i) == 0)
        
        if (nargin == 3) && ~isempty(start_point)
            continue
        else
            start_point = empty_track_start_point(1);
            start_point(1,sp_x0_idx) = x(3,i);
            start_point(1,sp_y0_idx) = x(4,i);
            start_point(1,sp_phi0_idx) = atan2d(x(4,i)-x(2,i),x(3,i)-x(1,i));
            i = i+1;
        end % if        
    
    % Unkown case
    else
        
        error('X2ContRailwayMap: Unhandled case!');
        
    end % if
    
end % while i

% Set starting point ______________________________________________________
if isempty(start_point)
    start_point = empty_track_start_point(1);
    start_point(1,sp_x0_idx) = x(1,1);
    start_point(1,sp_y0_idx) = x(2,1);
    start_point(1,sp_phi0_idx) = atan2d(x(4,1)-x(2,1),x(3,1)-x(1,1));
end % if

% if any(x_descriptor == 11)
%     start_point = empty_track_start_point(1);
% end % if

% Init variables __________________________________________________________
x = x(:,(x_descriptor~=0 & x_descriptor~=99));
track_maps = repmat({empty_track_map(nan)},size(track_elements,1),1);
track_maps_cov = nan(3,3,size(track_elements,1));
track_start_points = repmat(empty_track_start_point(nan),size(track_elements,1),1);
track_start_points_cov = nan(3,3,size(track_elements,1));

%% Calculations II/II

% Write to continuous railway map _________________________________________

x_column_index = 1;
for i = 1:size(track_elements,1)
    
    % Init ________________________________________________________________
    
    % Set railway-map start point
    if (i == start_point(:,sp_id_idx))
        track_start_points(i,:) = start_point;
    elseif isnan(track_start_points(i,sp_id_idx))
        track_start_points(i,:) = empty_track_start_point(i);
    end % if
    
    % Append topology
    topology = [ ... 
                 zeros(i-1,1),    eye(i-1); ...
                   zeros(1,1),zeros(1,i-1)  ...
               ]; 
    
	% Create continous railway-map ________________________________________
    
    track_element_i = track_elements(i);
    track_map_i = empty_track_map(nan);
    switch track_element_i
        
        case 1 % straight
            
           x_0 = x(1,x_column_index);
           y_0 = x(2,x_column_index);
           x_end = x(3,x_column_index);
           y_end =x(4,x_column_index);
           straight_length = norm([x_0;y_0]-[x_end;y_end]);

           track_map_i(1,tm_id_idx) = i;
           track_map_i(1,tm_element_idx) = 1;
           track_map_i(1,tm_r0_idx) = 0;
           track_map_i(1,tm_rend_idx) = 0;
           track_map_i(1,tm_length_idx) = straight_length;
           track_map_i(1,tm_vmax_idx) = calcSpeedLimit(v_max,0);           
                      
           x_column_index = x_column_index + 1;
           
        case 11 % transition straight 
            
           x_0 = x(3,x_column_index-1);
           y_0 = x(4,x_column_index-1);
           x_end = x(1,x_column_index+1);
           y_end = x(2,x_column_index+1);
           x_end_suc_track = x(3,x_column_index+1);
           y_end_suc_track = x(4,x_column_index+1);
           straight_length = norm([x_0;y_0]-[x_end;y_end]);
           phi_0 = atan2d(y_end-y_0,x_end-x_0);
           phi_0_suc_track = atan2d(y_end_suc_track-y_end,x_end_suc_track-x_end);

           track_map_i(1,tm_id_idx) = i;
           track_map_i(1,tm_element_idx) = 11;
           track_map_i(1,tm_r0_idx) = 0;
           track_map_i(1,tm_rend_idx) = 0;
           track_map_i(1,tm_length_idx) = straight_length;
           track_map_i(1,tm_vmax_idx) = calcSpeedLimit(v_max,0);
           
           track_start_points(i,:) = empty_track_start_point(i);
%            track_start_points(i,sp_x0_idx) = x_0;
%            track_start_points(i,sp_y0_idx) = y_0;
           track_start_points(i,sp_phi0_idx) = phi_0;
           
           track_start_points(i+1,:) = empty_track_start_point(i+1);
%            track_start_points(i+1,sp_x0_idx) = x_end;
%            track_start_points(i+1,sp_y0_idx) = y_end;
           track_start_points(i+1,sp_phi0_idx) = phi_0_suc_track;
                                
           x_column_index = x_column_index + 1;
           
        case 2 % clothoid (from straight to arc)
            
           clothoid_length = x(1,x_column_index);
           clothoid_radius = x(2,x_column_index);
            
           track_map_i(1,tm_id_idx) = i;
           track_map_i(1,tm_element_idx) = 2;
           track_map_i(1,tm_r0_idx) = 0; 
           track_map_i(1,tm_rend_idx) = clothoid_radius;
           track_map_i(1,tm_length_idx) = clothoid_length;
           track_map_i(1,tm_vmax_idx) = calcSpeedLimit(v_max,clothoid_radius);
           
        case 3 % circular arc
            
           arc_radius = x(2,x_column_index);
           arc_length = x(3,x_column_index);           
            
           track_map_i(1,tm_id_idx) = i;
           track_map_i(1,tm_element_idx) = 3;
           track_map_i(1,tm_r0_idx) = arc_radius; 
           track_map_i(1,tm_rend_idx) = arc_radius;
           track_map_i(1,tm_length_idx) = arc_length;
           track_map_i(1,tm_vmax_idx) = calcSpeedLimit(v_max,arc_radius);
           
           if x(4,x_column_index) == 0
               x_column_index = x_column_index + 1;
           end % if           

        case 4 % clothoid (from arc to straight)
            
           clothoid_radius = x(2,x_column_index);
           clothoid_length = x(4,x_column_index);
            
           track_map_i(1,tm_id_idx) = i;
           track_map_i(1,tm_element_idx) = 4;
           track_map_i(1,tm_r0_idx) = clothoid_radius;
           track_map_i(1,tm_rend_idx) = 0;
           track_map_i(1,tm_length_idx) = clothoid_length;
           track_map_i(1,tm_vmax_idx) = calcSpeedLimit(v_max,clothoid_radius);
           
           x_column_index = x_column_index + 1;
           
        case 5 % clothoid (from arc to arc)
            
           clothoid_start_radius = x(2,x_column_index);
            
           x_column_index = x_column_index + 1;
            
           clothoid_length = x(1,x_column_index);           
           clothoid_end_radius = x(2,x_column_index);
                       
           track_map_i(1,tm_id_idx) = i;
           track_map_i(1,tm_element_idx) = 5;
           track_map_i(1,tm_r0_idx) = clothoid_start_radius;
           track_map_i(1,tm_rend_idx) = clothoid_end_radius;
           track_map_i(1,tm_length_idx) = clothoid_length;
           track_map_i(1,tm_vmax_idx) = calcSpeedLimit(v_max,max(abs(clothoid_start_radius),abs(clothoid_end_radius)));
           
    end % switch
    
	track_maps{i} = track_map_i;
    
end % for i

%% Output

if any(isnan(track_start_points(:,sp_x0_idx)))
        
    num_sp = size(track_start_points,1);
    track_start_points_old = track_start_points;
    
    [~,~,track_start_points,~] = ... 
        calcMatTrackStartPoints([],topology,track_start_points,track_maps,track_start_points_cov,track_maps_cov);
    
    % Respect changed starting points because of transition straights _____
    
    if any(x_descriptor==11)
        track_maps_temp = cell2mat(track_maps);  
        trans_straight_idx = find(track_maps_temp(:,tm_element_idx)==11);
        for i = 1:length(trans_straight_idx)
            track_idx_i = trans_straight_idx(i);        
            track_id_i = track_maps_temp(track_idx_i,tm_id_idx);
            sp_idx_i = find(ismember(track_start_points(:,sp_id_idx),track_id_i));

            track_length_i = track_maps_temp(track_idx_i,tm_length_idx);
            track_p_0_i = [ track_start_points(sp_idx_i,sp_x0_idx); ... 
                            track_start_points(sp_idx_i,sp_y0_idx)];
            track_phi_0_i = track_start_points_old(sp_idx_i,sp_phi0_idx);

            suc_track_sp = straightLineTrackElement(track_length_i,track_phi_0_i,track_p_0_i);
            suc_track_phi_0 = track_start_points_old(sp_idx_i+1,sp_phi0_idx);

            % Set updated start-points around transition straight
            track_start_points(sp_idx_i,sp_phi0_idx) = track_start_points_old(sp_idx_i,sp_phi0_idx);
            track_start_points(sp_idx_i+1,sp_x0_idx) = suc_track_sp(1);
            track_start_points(sp_idx_i+1,sp_y0_idx) = suc_track_sp(2);
            track_start_points(sp_idx_i+1,sp_phi0_idx) = suc_track_phi_0;

            % Update start-points
            if (sp_idx_i+2) <= num_sp
                track_start_points(sp_idx_i+2:end,2:end) = nan;

                [~,~,track_start_points,~] = ... 
                    calcMatTrackStartPoints([],topology,track_start_points,track_maps,track_start_points_cov,track_maps_cov);
            end % if
        end % for i
    end % if
    
end % if

track_maps = cell2mat(track_maps);  

end % function

%% Helper Functions

function x_restored = restoreX(x,x_descriptor)

% Respect transition straights / fit x to length of x_descriptor
x_restored = zeros(4,sum(x_descriptor~=99));
x_idx = 1;
for i = 1:length(x_descriptor)
    
    if x_descriptor(i) == 11
        continue
    end % if
    
    if (x_descriptor(i) == 99)
        continue
    end % if
    
    x_restored(:,i) = x(:,x_idx);
    x_idx = x_idx + 1;
    
end % for i

end

