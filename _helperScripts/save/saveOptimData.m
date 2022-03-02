% 
% Save output data from 'localization.slx'
%   
%   Author: Hanno Winter
%   Date: 10-Apr-2021; Last revision: 18-Jul-2021

%% Init

if ~exist('optimout_data_saved','var') || optimout_data_saved == false
    optimout_data_saved = false;
    fprintf('Save optimization output data...');
elseif optimout_data_saved
    fprintf('Already saved optimization output data!\n');
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
    specifier_str = 'optimout_';
else
    specifier_str = 'optimout_atf_'; % 'atf' like artificial
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

%% Save

regexp_str = '.*optim.*';
save_vars_info = whos('-regexp',regexp_str);

if any([save_vars_info.bytes] > 2^31)
    save(save_name,'-regexp',regexp_str,'-v7.3')
else
    save(save_name,'-regexp',regexp_str)
end % if

%% Finish

optimout_data_saved = true;
fprintf('done!\n');
