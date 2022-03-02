% 
% Save reference track-map related to the IMM data (only relevant for data set Chemnitz)
%   
%   Author: Hanno Winter
%   Date: 13-Apr-2021; Last revision: 09-May-2021

%% Init

if ~exist('map_data_saved','var') || map_data_saved == false
    map_data_saved = false;
    fprintf('Save reference track-map data...');
elseif map_data_saved
    fprintf('Already saved reference track-map data!\n');
    return  
end % if

%% Settings

% Set Output file name
%   Current pattern: 'simout_<data set>_<start time>-<end time>-<date>'

if ~bdIsLoaded(sim_model_name)
    %load_system(sim_model_name);
    sim_stop_time_str = 'yyyy-MM-dd HH:mm:ss';
else
    sim_duration = str2double(get_param(sim_model_name, 'StopTime'));
    sim_stop_time_str = datetime(start_time,'InputFormat','yyyy-MM-dd HH:mm:ss') + seconds(sim_duration);
    sim_stop_time_str = datestr(sim_stop_time_str,'yyyy-mm-dd HH:MM:SS');
end % if

if ~insert_gps_outages
    specifier_str = 'map_data_';
else
    specifier_str = 'map_data_atf_'; % 'atf' like artificial
end % if

save_name = [ ... 
              specifier_str, ... 
              input_data_selector, ... 
              '_', ... 
              num2str(start_time([12,13,15,16,18,19])), ...
              '-', ...
              num2str(sim_stop_time_str([12,13,15,16,18,19])), ...
              '-', ...
              num2str(start_time([9,10,6,7,3,4])) ... 
            ];

%% Calc

save( ... 
      save_name, ... 
      'map_data_prepared', ...
      'imm_ref_track_map', ... 
      'plot_ref_map', ... 
      'ref_track_map_selector', ... 
      'ref_track_map', ... 
      'exclude_stillstand_cdf', ...
      'v_cdf_min', ... 
      'harmonize_with_gps' ... 
    );

%% Finish

map_data_saved = true;
fprintf('done!\n');

