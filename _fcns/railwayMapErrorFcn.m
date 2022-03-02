function F = railwayMapErrorFcn(x,x_descriptor,z,R,z_data_selector,config)
% F = railwayMapErrorFcn(x,x_descriptor,z,R,z_data_selector)
% 

%% Settings

straight_sloppiness_threshold = config.straight_sloppiness.threshold; % in m (% Start punishing length of straights when it is getting longer than this value (in m), compared to the initial length)
straight_sloppiness_range = config.straight_sloppiness.range; % in m (% Ramp length around straight sloppiness threshold to smooth punishing too long straights)
min_z_per_straight = config.min_z_per_element.straight; % Dismiss straights with less measurements
min_z_per_arc = config.min_z_per_element.arc; % Dismiss circular-arcs with less measurements

%% Initialization

if all(isnan(x(:)))
    error('railwayErrorFcn: ''x'' got messed up!')
end % if

% Set some static numbers and codes
straight_code = 1;
arc_code = 3;

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

% Reshape 
x = reshape(x,4,length(x(:))/4);

% Check for negative length values
x = ensurePositiveLength(x,x_descriptor);

% % Limit radius
% x = ensureRadiusLimits(x,x_descriptor);

% Convert to continous railway-map
[~,track_maps,track_start_points] = x2contMatRailwayMap(x,x_descriptor,[]);

%% Calculations 

% Calculate error vectors and information matrices ________________________

[r_straight,r_straight_lengths,omega_r_straight] = calcErrorToStraightGeometry(z,R,z_data_selector,track_maps,track_start_points);
r_straight_empty_selector = cellfun(@(cell) isempty(cell) | (size(cell,2) < min_z_per_straight),r_straight);
r_straight = r_straight(~r_straight_empty_selector);
omega_r_straight = omega_r_straight(~r_straight_empty_selector);

[r_arc,r_arc_lengths,omega_r_arc] = calcErrorToCircularGeometry(z,R,z_data_selector,track_maps,track_start_points);
r_arc_empty_selector = cellfun(@(cell) isempty(cell) | (size(cell,2) < min_z_per_arc),r_arc);
r_arc = r_arc(~r_arc_empty_selector);
omega_r_arc = omega_r_arc(~r_arc_empty_selector);

% Create F_i ______________________________________________________________

F_r_straight_length = cell2mat(r_straight_lengths(~r_straight_empty_selector));
F_r_straight_length = reshape(F_r_straight_length,[2,length(F_r_straight_length)/2]);
F_r_straight_length = sum(F_r_straight_length,1);
F_r_straight_length = F_r_straight_length(:);

F_r_straight = zeros(sum(cellfun(@(cell) size(cell,2),r_straight)),1);
for i = 1:size(r_straight,1)    
    
    F_r_straight_i = zeros(1,size(r_straight{i},2));
    for j = 1:size(r_straight{i},2)
        F_r_straight_i(:,j) = sqrt(r_straight{i}(:,j)'*omega_r_straight{i}(:,:,j)*r_straight{i}(:,j));
%         F_r_straight_i(:,j) = sqrt(r_straight{i}(:,j)'*r_straight{i}(:,j));
        if any(isnan(F_r_straight_i(:))) || any(isinf(F_r_straight_i(:)))
            F_r_straight_i(isnan(F_r_straight_i)) = 0;
            F_r_straight_i(isinf(F_r_straight_i)) = 0;
        end % if
    end % for j
    
    if i > 1
        start_index = sum(cellfun(@(cell) size(cell,2),r_straight(1:i-1))) + 1;
    else
        start_index = 1;
    end % if
    end_index = start_index -1 + size(r_straight{i},2);  
    F_r_straight(start_index:end_index) = F_r_straight_i;
    
end % for k

F_r_arc_length = cell2mat(r_arc_lengths(~r_arc_empty_selector));

F_r_arc = zeros(sum(cellfun(@(cell) size(cell,2),r_arc)),1);
for i = 1:size(r_arc,1)
    
    F_r_arc_i = zeros(1,size(r_arc{i},2));
    for j = 1:size(r_arc{i},2)
        F_r_arc_i(:,j) = sqrt(r_arc{i}(:,j)'*omega_r_arc{i}(:,:,j)*r_arc{i}(:,j));
%         F_r_arc_i(:,j) = sqrt(r_arc{i}(:,j)'*r_arc{i}(:,j));
        if any(isnan(F_r_arc_i(:))) || any(isinf(F_r_arc_i(:)))
            F_r_arc_i(isnan(F_r_arc_i)) = 0;
            F_r_arc_i(isinf(F_r_arc_i)) = 0;
        end % if
    end % for j  
    
    if i > 1
        start_index = sum(cellfun(@(cell) size(cell,2),r_arc(1:i-1))) + 1;
    else
        start_index = 1;
    end % if
    end_index = start_index -1 + size(r_arc{i},2);  
    F_r_arc(start_index:end_index) = F_r_arc_i;
    
end % for k

% Calculate error from track-start and end to first/last measurment point _

% Get first measurement
first_track_id = track_maps(1,tm_id_idx);
first_track_element = track_maps(1,tm_element_idx);
first_z_idx = find(z_data_selector(z_data_selector(:,1)==first_track_id,2:end),1,'first');
z_first = z(:,first_z_idx);

% Get first track position
tsp_idx = (track_start_points(:,sp_id_idx)==first_track_id);
x_0 = track_start_points(tsp_idx,sp_x0_idx);
y_0 = track_start_points(tsp_idx,sp_y0_idx);
p_first = [x_0;y_0];

% Get last measurment
last_track_id = track_maps(end,tm_id_idx);
last_track_element = track_maps(end,tm_element_idx);
last_z_idx = find(z_data_selector(z_data_selector(:,1)==last_track_id,2:end),1,'last');
z_last = z(:,last_z_idx);

% Get last track position
track_length = track_maps(end,tm_length_idx);
track_radius = track_maps(end,tm_rend_idx);
tsp_idx = (track_start_points(:,sp_id_idx)==last_track_id);
x_0 = track_start_points(tsp_idx,sp_x0_idx);
y_0 = track_start_points(tsp_idx,sp_y0_idx);
track_p_0 = [x_0;y_0];
phi_0 = track_start_points(tsp_idx,sp_phi0_idx);
switch last_track_element
    case straight_code
        [p_last,~,~,~] = straightLineTrackElement(track_length,phi_0,track_p_0);
    case arc_code
        [p_last,~,~,~] = constantArcTrackElement(track_length,track_radius,phi_0,track_p_0);
    otherwise
        p_last = [];
end % if

% Calculate error between first z and first track-end
if ~isempty(p_first) && ~isempty(z_first)
    F_r_first_point = norm(z_first-p_first);
%     if F_r_first_point < 100
%         F_r_first_point = 0;
%     end % if
else
    F_r_first_point = [];
end % if

% Calculate error between last z and last track-end
if ~isempty(p_last) && ~isempty(z_last)
    F_r_end_point = sqrt(norm(z_last-p_last));
%     if F_r_end_point < 100
%         F_r_end_point = 0;
%     end % if
else
    F_r_end_point = [];
end % if

%% Adjust

% Normalize _______________________________________________________________

% w_straight_desired = 0.5;
% w_arc_desired = 0.5;
% 
% n_straight = length(F_r_straight);
% n_arc = length(F_r_arc);
% n_all = n_straight + n_arc;
% 
% w_straight = w_straight_desired*n_all/n_straight;
% w_arc = w_arc_desired*n_all/n_arc;
% 
% if isinf(w_straight)
%     w_arc = 1;
% end % if
% 
% if isinf(w_arc)
%     w_straight = 1;
% end % if

% Allow some sloppiness in length of straights ____________________________

ramp = @(x) min(max( 1/straight_sloppiness_range*(x-(straight_sloppiness_threshold-straight_sloppiness_range/2)) , 0) , 1);
F_r_straight_length = ramp(F_r_straight_length).*F_r_straight_length;
F_r_straight_length = sqrt(F_r_straight_length);

% if any(F_r_straight_length)
%     test = 1;
% end % if

%% Output 

% F = [ ...
%       w_straight * F_r_straight; ... 
%            w_arc *      F_r_arc; ... 
%                   F_r_end_point  ... 
%     ];

% F = [ F_r_straight; F_r_arc; F_r_end_point ];
% F = [ F_r_straight; F_r_arc; F_r_end_point ];
F = [ F_r_straight; F_r_straight_length; F_r_arc; F_r_end_point ];
% F = [ F_r_straight; F_r_arc ];

end
