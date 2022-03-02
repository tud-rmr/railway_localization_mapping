%
%   Data available at: http://dx.doi.org/10.25534/tudatalib-360
%
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 02-Dec-2020
% 

%% Init

loadBsFilenamePatterns

%% Checks

% if ~exist('import_BS_map_data_executed','var')
%     import_BS_map_data_executed = false;
% elseif import_BS_map_data_executed
%     fprintf('Import map data script already executed!\n');
%     return  
% end % if

%% Map Data
%
%   Data available at: http://dx.doi.org/10.25534/tudatalib-360
%
                             
map_xy_track_map_data_paths = { ... 
                                fullfile('Dataset_BS_20190222','_maps','map.csv') ... 
                              };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map: (X,Y)-Track-Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[map_xy_track_map, ~, load_flag] = importFromFile(map_xy_track_map_data_paths,'ref_track_map_Bs','DataType','map');
if load_flag == 1 % save to .mat for faster access in the future
    writeToMatFile(map_xy_track_map,'ref_track_map_Bs','ref_track_map');
end % for i

%% Finish script

import_BS_map_data_executed = true;
