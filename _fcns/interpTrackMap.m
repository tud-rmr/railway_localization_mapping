function track_map_new = interpTrackMap(d,d_ref,track_map)
% track_map_new = interpTrackMap(d,d_ref,track_map)
% 
%   Author: Hanno Winter
%   Date: 11-Apr-2020; Last revision: 11-Apr-2020

%%

track_map_new = table();

variable_names = track_map.Properties.VariableNames;
for i = 1:length(variable_names)
    name_i = variable_names{i};
    value_i = track_map.(name_i);
    
    track_map_new.(name_i) = interp1(d_ref(:),value_i(:),d(:),'pchip');
end % if

end

