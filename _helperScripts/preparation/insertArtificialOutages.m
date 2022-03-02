% 
% Insert artificial GNSS outages
%   
%   Author: Hanno Winter
%   Date: 18-Jul-2021; Last revision: 18-Jul-2021

%% Settings

% Load settings
configSimulation

%% Insert artificial GNSS outages

if insert_gps_outages
    
    % Init ________________________________________________________________
    
    % Make sure to always start from the original validity data
    gnss_data_validity = gnss_data_validity_original; 
    
    % Random number generator seed
    rng_seed = sum(int8(input_data_selector)) + ... 
           round(target_gps_availability*100) + ... 
           round(target_gps_outage_length(1)) + ... 
           round(target_gps_outage_length(2));
    
    gnss_data_valid_indices = find(gnss_data_validity);
    
    % If simulink modal is already loaded, adapt calculations to simulation time
    if bdIsLoaded(sim_model_name)
        sim_stop_time = str2double(get_param(sim_model_name,'StopTime'));
        
        max_time_idx = find(time_gnss <= sim_stop_time,1,'last');
        % gnss_valid_max_index = find(gnss_data_valid_indices<=max_time_idx,1,'last');
    else
        max_time_idx = length(time_gnss);
        % gnss_valid_max_index = length(gnss_data_valid_indices);
    end % if  
    
    % Create temporary validity data
    gnss_data_validity_temp = gnss_data_validity(1:max_time_idx);
    gnss_data_valid_indices_temp = find(gnss_data_validity_temp);
    
    % Get quantity of random numbers needed
    n_gps_outages = sum(~gnss_data_validity_temp);
    n_gps_full = length(gnss_data_validity_temp);
    n_artificial_outages = (1-target_gps_availability)*n_gps_full-n_gps_outages;
    n_artificial_outages = round(n_artificial_outages / mean(target_gps_outage_length));
        
    % Determine outage lengths and times
    rng(rng_seed,'twister');    
    outage_length = randi(target_gps_outage_length,n_artificial_outages,1);
    
    rng(rng_seed,'twister');
    outage_times_indices_temp = randi(length(gnss_data_valid_indices_temp),n_artificial_outages,1);
    outage_times_indices = gnss_data_valid_indices_temp(outage_times_indices_temp);
    
    % Adapt valid GNSS data selector
    for i = 1:length(outage_times_indices)
        outage_start_index = outage_times_indices(i);
        outage_start_index = max(1,outage_start_index);
        outage_end_index = outage_start_index + ceil(outage_length(i)/natural_gnss_sample_time) - 1;
        outage_end_index = min(length(gnss_data_validity_temp),outage_end_index);
        
        outage_mask = true(size(gnss_data_validity));
        outage_mask(outage_start_index:outage_end_index) = false;
                
        gnss_data_validity = gnss_data_validity & outage_mask;
    end % for i

    % Update validity data ________________________________________________
    gnss_sim_data_validity = timeseries(gnss_data_validity,time_gnss');
    
    % Report ______________________________________________________________
    new_gnss_availability = sum(gnss_data_validity(1:max_time_idx))/length(gnss_data_validity_temp);
    old_gnss_availability = sum(gnss_data_validity_original(1:max_time_idx))/length(gnss_data_validity_temp);
    fprintf('\nInserted artificial GNSS outages. GNSS availability is now: %.2f (instead of %.2f)\n',new_gnss_availability,old_gnss_availability);
    
end % if
