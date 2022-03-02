% 
%   Data available at: http://dx.doi.org/10.25534/tudatalib-166.2
%
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 18-Apr-2021
% 

%% Init

loadCFilenamePatterns

%% Checks

% if ~exist('import_C_map_data_executed','var')
%     import_C_map_data_executed = false;
% elseif import_C_map_data_executed
%     fprintf('Import map data script already executed!\n');
%     return  
% end % if

%% Map Data
% 
%   Data available at: http://dx.doi.org/10.25534/tudatalib-166.2
%
                             
ref_track_map_data_paths = { ... 
                                fullfile('Dataset_C_20181024','_maps','ref_track_map.csv') ... 
                              };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map: (X,Y)-Track-Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ref_track_map, ~, load_flag] = importFromFile(ref_track_map_data_paths,'ref_track_map_C','DataType','map');
if load_flag == 1 % save to .mat for faster access in the future
    writeToMatFile(ref_track_map,'ref_track_map_C','ref_track_map');
end % for i

%% Finish script

import_C_map_data_executed = true;
