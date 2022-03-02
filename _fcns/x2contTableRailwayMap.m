function railway_map = x2contTableRailwayMap(x,x_descriptor,railway_map_start_point)
% railway_map = x2contTableRailwayMap(x,x_descriptor,railway_map_start_point)
% 

%% Init

if (nargin < 3)
    railway_map_start_point = [];
elseif istable(railway_map_start_point)
    railway_map_start_point = tableTrackStartPoints2matTrackStartPoints(railway_map_start_point);
end % if

%% Calculations

[topology,track_maps,track_start_points] = x2contMatRailwayMap(x,x_descriptor,railway_map_start_point);

railway_map.topology = topology;
railway_map.track_start_points = matTrackStartPoints2tableTrackStartPoints(track_start_points);
railway_map.track_maps = matTrackMap2tableTrackMap(track_maps);

end % function

