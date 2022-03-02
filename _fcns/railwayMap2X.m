function [x,x_descriptor] = railwayMap2X(railway_map)
% [x,x_descriptor] = railwayMap2X(railway_map)
% 

x = zeros(4,1);
x_descriptor = [];
x_column_index = 1;
for i = 1:size(railway_map.track_maps,1)
    
    % Identify current track-element, as well as the surrounding ones _____
    
    h = i-1;
    if h > 0
        track_element_h = railway_map.track_maps(h,:).track_element;
    else
        track_element_h = 0;
    end % if 
    
    track_element_i = railway_map.track_maps(i,:).track_element;
    
    j = i+1;
    if j <= size(railway_map.track_maps,1)
         track_element_j = railway_map.track_maps(j,:).track_element;
    else
        track_element_j = 0;
    end % if   
        
    % Build vector railway-map depending on current track-element
    switch track_element_i
        
        case 1 % straight
            
            % Consider orientation of straight
            if abs(railway_map.track_maps(i,:).length) > 1e-2
                x_0 = railway_map.track_start_points(i,:).x_0;
                y_0 = railway_map.track_start_points(i,:).y_0;
                track_end_point = calcTrackEndPoints(railway_map,railway_map.track_start_points(i,:).ID);
                x_end = track_end_point.x_end;
                y_end = track_end_point.y_end;
            else
                phi_0 =  railway_map.track_start_points(i,:).phi_0;
                
                x_end = railway_map.track_start_points(i,:).x_0;
                y_end = railway_map.track_start_points(i,:).y_0;
                x_0 = x_end - cosd(phi_0);
                y_0 = y_end - sind(phi_0);
            end % if
            
            x(1,x_column_index) = x_0;
            x(2,x_column_index) = y_0;
            x(3,x_column_index) = x_end;
            x(4,x_column_index) = y_end;
            
            x_descriptor(1,x_column_index) = 1;
            
            x_column_index = x_column_index + 1;
        
        case 11 % transition straight
            
            x_descriptor(1,x_column_index) = 11;
            x_column_index = x_column_index + 1;
            
        case 2 % clothoid (from straight to arc)
            
            x(1,x_column_index) = railway_map.track_maps(i,:).length;
            %x(2:4,x_column_index) = 0;
            
        case 3 % circular arc
            
            x(2,x_column_index) = railway_map.track_maps(i,:).r_end;
            x(3,x_column_index) = railway_map.track_maps(i,:).length;
            % x(4,x_column_index) = 0;
            
            circle_pos_str = [num2str(track_element_h),num2str(track_element_i),num2str(track_element_j)];
            
            % Consider surroundings
            switch circle_pos_str
                case {'131'} % surrounded by straights
                    x_descriptor(1,x_column_index) = 131;
                    x_column_index = x_column_index + 1;
                case {'134','135','130'} % previous element is a straight
                    x_descriptor(1,x_column_index) = 13;
                case {'231','531','031'} % succeeding element is a straight
                    x_descriptor(1,x_column_index) = 31;
                    x_column_index = x_column_index + 1;
                otherwise % arc surrounded by clothoids
                    x_descriptor(1,x_column_index) = 3;
            end % switch
            
        case 4 % clothoid (from arc to straight)
            
            x(4,x_column_index) = railway_map.track_maps(i,:).length;
            
            x_column_index = x_column_index + 1;
            
        case 5 % clothoid (from arc to arc)
            
            x(4,x_column_index) = railway_map.track_maps(i,:).length;
            
            x_column_index = x_column_index + 1;
            
            x(1,x_column_index) = railway_map.track_maps(i,:).length;
            %x(2:4,x_column_index) = 0;
            
    end % switch
    
end % for i

% Remove none necessary entries for transion straights
zero_column_selector = (sum(x,1) == 0);
x = x(:,~zero_column_selector);

% Set start-point explicitly, when first element is an arc
if contains(num2str(x_descriptor(1)),'3')
    x_descriptor = [0, x_descriptor];
    
    x_0_start = railway_map.track_start_points(1,:).x_0;
    y_0_start = railway_map.track_start_points(1,:).y_0;
    phi_0_start = railway_map.track_start_points(1,:).phi_0;
    
    x_start(1,1) = x_0_start - cosd(phi_0_start);
    x_start(2,1) = y_0_start - sind(phi_0_start);
    x_start(3,1) = x_0_start;
    x_start(4,1) = y_0_start;
    
    x = [x_start,x];
end % if

% Set end-point explicitly, when last element is an arc
if contains(num2str(x_descriptor(end)),'3')
    x_descriptor = [x_descriptor,99];
    
%     [track_end_point,~] = calcTrackEndPoints(railway_map,railway_map.track_start_points(end,:).ID);
%         
%     x_end(1,1) = track_end_point.x_end - cosd(track_end_point.phi_end);
%     x_end(2,1) = track_end_point.y_end - sind(track_end_point.phi_end);
%     x_end(3,1) = track_end_point.x_end;
%     x_end(4,1) = track_end_point.y_end;
%     
%     x = [x,x_end];
end % if

end % function

