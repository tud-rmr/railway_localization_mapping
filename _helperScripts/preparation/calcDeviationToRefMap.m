% 
% Calculate deviation to reference-map
%   
%   Author: Hanno Winter
%   Date: 25-Apr-2021; Last revision: 05-Jul-2021

%% Settings

v_cdf_min = v_min_tgc_speed-v_min_tgc_hysteresis;

%% Init

if ~ismember(input_data_selector,{'C','BS'})
    return
end % if

if exist('run_optimization_executed','var') && run_optimization_executed
    prepareOptimOutData
end % if

fprintf('Calculating deviations to reference-map:\n')

ref_map_abs_position_utm_orig = [ref_track_map.UtmEast_m(ref_track_map_selector)';ref_track_map.UtmNorth_m(ref_track_map_selector)']-[p_0_utm(1);p_0_utm(2)];
ref_map_abs_position_utm = ref_map_abs_position_utm_orig;

% % CDF time vector _________________________________________________________
% 
% t_cdf = simulink_time(1):ts_cdf:simulink_time(end);

% Exclude parts with slow speed ___________________________________________

v_min_selector = abs(sim_imm.VelocityVehicle_ms) > v_cdf_min;
v_min_gps_selector = v_min_selector(simout_gnss_valid_flag.Data(:,1)==1);

% Exclude parts where vehicle was not on the reference track ______________

% Station Schlettau
exclusion_zone1_ll = [ 50.5535270, 50.556035; ...
                       12.9509635, 12.957117];
[exclusion_zone1_utm_x,exclusion_zone1_utm_y] = ll2utm(exclusion_zone1_ll(1,:)',exclusion_zone1_ll(2,:)');
exclusion_zone1_utm = [exclusion_zone1_utm_x(:)';exclusion_zone1_utm_y(:)']-[p_0_utm(1);p_0_utm(2)];

exclusion_zone1_selector = isInArea(ref_map_abs_position_utm,exclusion_zone1_utm(1,1),exclusion_zone1_utm(1,2),exclusion_zone1_utm(2,1),exclusion_zone1_utm(2,2));
exclusion_zones_selector = exclusion_zone1_selector;

% ref_map_abs_position_utm = ref_map_abs_position_utm(:,~exclusion_zones_selector);
ref_map_abs_position_utm(:,exclusion_zones_selector) = nan;

%% GNSS

fprintf('\tGNSS data...')

% Availability ____________________________________________________________

expected_gnss_time = simulink_time(1):natural_gnss_sample_time:simulink_time(end);
num_expected_gnss = simulink_time(end)/natural_gnss_sample_time+1; % +1 to include first time step at t=0
num_valid_gnss = length(sim_gnss_filter_input.Time);
num_missing_gnss_z = num_expected_gnss-num_valid_gnss;
sim_gnss_availability = num_valid_gnss/num_expected_gnss; 

% Perpendicular deviation from ref-map ____________________________________

% Original data
p_test_orig = [sim_gnss_filter_input.UtmEast_m';sim_gnss_filter_input.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
time_orig = sim_gnss_filter_input.Time;
p_test = p_test_orig;

% Exclude parts where vehicle is not on reference track
exclusion_zone1_selector = isInArea(p_test,exclusion_zone1_utm(1,1),exclusion_zone1_utm(1,2),exclusion_zone1_utm(2,1),exclusion_zone1_utm(2,2));
exclusion_selector = exclusion_zone1_selector;
p_test(:,exclusion_selector) = nan;

% Further data selection
if exclude_stillstand_cdf    
    num_missing_gnss_z = ceil(sum(v_min_gps_selector)*(1-sim_gnss_availability));    
    p_test = p_test(:,v_min_gps_selector');
end % if
v_min_gps_indices = find(v_min_gps_selector);

d_perpendicular_e_pp = nan(size(time_orig));
[d_perpendicular_e_pp_temp,~,p_ref_e_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm_orig,p_test_orig);
d_perpendicular_e_pp(p_ref_e_pp_indices) = d_perpendicular_e_pp_temp;

% [d_perpendicular,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm(:,~isnan(ref_map_abs_position_utm(1,:))),p_test);
[d_perpendicular,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,p_test);
d_perpendicular = abs(d_perpendicular);
d_perpendicular_temp = [d_perpendicular,nan(1,num_missing_gnss_z)];
p_ref_pp = [p_ref_pp,nan(2,num_missing_gnss_z)];

% PP Deviation from map CDF _______________________________________________

e_pp_cdf_error = 0:0.01:ceil(max(d_perpendicular)); % ensure to include all data by including maximum existing error
e_pp_cdf = [e_pp_cdf_error(2:end); zeros(1,length(e_pp_cdf_error)-1)];
e_pp_cdf(2,:) = histcounts(d_perpendicular_temp,e_pp_cdf_error,'Normalization','cdf');

% Track-length deviation CDF ______________________________________________

time_cdf = sim_gnss_filter_input.Time(v_min_gps_indices(p_ref_pp_indices));
none_nan_data_selector = ~isnan(p_ref_pp(1,:));
d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp(:,none_nan_data_selector),1,2)].^2,1)));
p_test_temp = p_test(:,p_ref_pp_indices);
p_test_temp = [p_test_temp,nan(2,num_missing_gnss_z)];
d_test = cumsum(sqrt(sum([zeros(2,1),diff(p_test_temp(:,none_nan_data_selector),1,2)].^2,1)));
e_track_length_orig = (d_p_ref_pp-d_test);
e_track_length = abs(e_track_length_orig);

e_tl_cdf_error = floor(min(e_track_length)):0.01:ceil(max(e_track_length)); % ensure to include all data by including maximum existing error
e_tl_cdf = [e_tl_cdf_error(2:end); zeros(1,length(e_tl_cdf_error)-1)];
e_tl_cdf(2,:) = histcounts(e_track_length,e_tl_cdf_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

% sim_gnss_filter_input_error_data.e_time = time_orig;
sim_gnss_filter_input_error_data.cdf_time = time_cdf;
sim_gnss_filter_input_error_data.e_pp = d_perpendicular_e_pp;
sim_gnss_filter_input_error_data.e_pp_cdf.MaxError_m = e_pp_cdf(1,:)';
sim_gnss_filter_input_error_data.e_pp_cdf.Availability = e_pp_cdf(2,:)';
sim_gnss_filter_input_error_data.e_tl = e_track_length_orig;
sim_gnss_filter_input_error_data.e_tl_cdf.MaxError_m = e_tl_cdf(1,:)';
sim_gnss_filter_input_error_data.e_tl_cdf.Availability = e_tl_cdf(2,:)';

fprintf('done!\n')

%% EKF

fprintf('\tEKF data...')

% Perpendicular deviation from ref-map ____________________________________

% Original data
p_test_orig = [sim_ekf.UtmEast_m';sim_ekf.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
time_orig = sim_ekf.Time;
p_test = p_test_orig;
d_test = sim_ekf.DistanceVehicle_m';

% Exclude parts where vehicle is not on reference track
exclusion_zone1_selector = isInArea(p_test,exclusion_zone1_utm(1,1),exclusion_zone1_utm(1,2),exclusion_zone1_utm(2,1),exclusion_zone1_utm(2,2));
exclusion_selector = exclusion_zone1_selector;
p_test(:,exclusion_selector) = nan;
d_test(exclusion_selector) = nan;

% Further data selection
time_selector = ismember(sim_ekf.Time,expected_gnss_time);
if exclude_stillstand_cdf && harmonize_with_gps
    data_selector = time_selector & v_min_selector;
elseif exclude_stillstand_cdf
    data_selector = v_min_selector;
elseif harmonize_with_gps
    data_selector = time_selector;
else
    data_selector = true(size(v_min_selector));
end % if
data_indices = find(data_selector);
p_test = p_test(:,data_selector);
d_test = d_test(data_selector);

% Calculations
d_perpendicular_e_pp = nan(size(time_orig));
[d_perpendicular_e_pp_temp,~,p_ref_e_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm_orig,p_test_orig);
d_perpendicular_e_pp(p_ref_e_pp_indices) = d_perpendicular_e_pp_temp;

[d_perpendicular,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,p_test);
d_perpendicular = abs(d_perpendicular);

% PP Deviation from map CDF _______________________________________________

e_pp_cdf_error = 0:0.01:ceil(max(d_perpendicular)); % ensure to include all data by including maximum existing error
e_pp_cdf = [e_pp_cdf_error(2:end); zeros(1,length(e_pp_cdf_error)-1)];
e_pp_cdf(2,:) = histcounts(d_perpendicular,e_pp_cdf_error,'Normalization','cdf');

% Track-length deviation CDF ______________________________________________

time_cdf = sim_ekf.Time(data_indices(p_ref_pp_indices));
none_nan_data_selector = ~isnan(p_ref_pp(1,:));
d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp(:,none_nan_data_selector),1,2)].^2,1)));
p_test_temp = p_test(:,p_ref_pp_indices);
d_test = cumsum(sqrt(sum([zeros(2,1),diff(p_test_temp(:,none_nan_data_selector),1,2)].^2,1)));
% d_test = d_test(p_ref_pp_indices);
d_test = d_test - d_test(1);
e_track_length_orig = (d_p_ref_pp-d_test);
e_track_length = abs(e_track_length_orig);

e_tl_cdf_error = floor(min(e_track_length)):0.01:ceil(max(e_track_length)); % ensure to include all data by including maximum existing error
e_tl_cdf = [e_tl_cdf_error(2:end); zeros(1,length(e_tl_cdf_error)-1)];
e_tl_cdf(2,:) = histcounts(e_track_length,e_tl_cdf_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

% sim_ekf_error_data.e_time = time_orig;
sim_ekf_error_data.cdf_time = time_cdf;
sim_ekf_error_data.e_pp = d_perpendicular_e_pp;
sim_ekf_error_data.e_pp_cdf.MaxError_m = e_pp_cdf(1,:)';
sim_ekf_error_data.e_pp_cdf.Availability = e_pp_cdf(2,:)';
sim_ekf_error_data.e_tl = e_track_length_orig;
sim_ekf_error_data.e_tl_cdf.MaxError_m = e_tl_cdf(1,:)';
sim_ekf_error_data.e_tl_cdf.Availability = e_tl_cdf(2,:)';

fprintf('done!\n')

%% IMM

fprintf('\tIMM data...')

% Perpendicular deviation from ref-map ____________________________________

% Original data
p_test_orig = [sim_imm.UtmEast_m';sim_imm.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
time_orig = sim_imm.Time;
p_test = p_test_orig;
d_test = sim_imm.DistanceVehicle_m';

% Exclude parts where vehicle is not on reference track
exclusion_zone1_selector = isInArea(p_test,exclusion_zone1_utm(1,1),exclusion_zone1_utm(1,2),exclusion_zone1_utm(2,1),exclusion_zone1_utm(2,2));
exclusion_selector = exclusion_zone1_selector;
p_test(:,exclusion_selector) = nan;
d_test(exclusion_selector) = nan;

% Further data selection
time_selector = ismember(sim_ekf.Time,expected_gnss_time);
if exclude_stillstand_cdf && harmonize_with_gps
    data_selector = time_selector & v_min_selector;
elseif exclude_stillstand_cdf
    data_selector = v_min_selector;
elseif harmonize_with_gps
    data_selector = time_selector;
else
    data_selector = true(size(v_min_selector));
end % if
data_indices = find(data_selector);
p_test = p_test(:,data_selector);
d_test = d_test(data_selector);

% Calculations
d_perpendicular_e_pp = nan(size(time_orig));
[d_perpendicular_e_pp_temp,~,p_ref_e_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm_orig,p_test_orig);
d_perpendicular_e_pp(p_ref_e_pp_indices) = d_perpendicular_e_pp_temp;

[d_perpendicular,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,p_test);
d_perpendicular = abs(d_perpendicular);

% PP Deviation from map CDF _______________________________________________

e_pp_cdf_error = 0:0.01:ceil(max(d_perpendicular)); % ensure to include all data by including maximum existing error
e_pp_cdf = [e_pp_cdf_error(2:end); zeros(1,length(e_pp_cdf_error)-1)];
e_pp_cdf(2,:) = histcounts(d_perpendicular,e_pp_cdf_error,'Normalization','cdf');

% Track-length deviation CDF ______________________________________________

time_cdf = sim_imm.Time(data_indices(p_ref_pp_indices));
none_nan_data_selector = ~isnan(p_ref_pp(1,:));
d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp(:,none_nan_data_selector),1,2)].^2,1)));
p_test_temp = p_test(:,p_ref_pp_indices);
d_test = cumsum(sqrt(sum([zeros(2,1),diff(p_test_temp(:,none_nan_data_selector),1,2)].^2,1)));
% d_test = d_test(p_ref_pp_indices);
d_test = d_test - d_test(1);
e_track_length_orig = (d_p_ref_pp-d_test);
e_track_length = abs(e_track_length_orig);

e_tl_cdf_error = floor(min(e_track_length)):0.01:ceil(max(e_track_length)); % ensure to include all data by including maximum existing error
e_tl_cdf = [e_tl_cdf_error(2:end); zeros(1,length(e_tl_cdf_error)-1)];
e_tl_cdf(2,:) = histcounts(e_track_length,e_tl_cdf_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

% sim_imm_error_data.e_time = time_orig;
sim_imm_error_data.cdf_time = time_cdf;
sim_imm_error_data.e_pp = d_perpendicular_e_pp;
sim_imm_error_data.e_pp_cdf.MaxError_m = e_pp_cdf(1,:)';
sim_imm_error_data.e_pp_cdf.Availability = e_pp_cdf(2,:)';
sim_imm_error_data.e_tl = e_track_length_orig;
sim_imm_error_data.e_tl_cdf.MaxError_m = e_tl_cdf(1,:)';
sim_imm_error_data.e_tl_cdf.Availability = e_tl_cdf(2,:)';

fprintf('done!\n')

%% TGC

fprintf('\tTGC data...')

% Perpendicular deviation from ref-map ____________________________________

% Original data
p_test_orig = [sim_tgc.UtmEast_m';sim_tgc.UtmNorth_m']-[p_0_utm(1);p_0_utm(2)];
time_orig = sim_tgc.Time;
p_test = p_test_orig;
d_test = sim_tgc.DistanceVehicle_m';

% Exclude parts where vehicle is not on reference track
exclusion_zone1_selector = isInArea(p_test,exclusion_zone1_utm(1,1),exclusion_zone1_utm(1,2),exclusion_zone1_utm(2,1),exclusion_zone1_utm(2,2));
exclusion_selector = exclusion_zone1_selector;
p_test(:,exclusion_selector) = nan;
d_test(exclusion_selector) = nan;

% Further data selection
time_selector = ismember(sim_ekf.Time,expected_gnss_time);
if exclude_stillstand_cdf && harmonize_with_gps
    data_selector = time_selector & v_min_selector;
elseif exclude_stillstand_cdf
    data_selector = v_min_selector;
elseif harmonize_with_gps
    data_selector = time_selector;
else
    data_selector = true(size(v_min_selector));
end % if
data_indices = find(data_selector);
p_test = p_test(:,data_selector);
d_test = d_test(data_selector);

p_test_optim_map_correspondence = p_test;

% Calculations
d_perpendicular_e_pp = nan(size(time_orig));
[d_perpendicular_e_pp_temp,~,p_ref_e_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm_orig,p_test_orig);
d_perpendicular_e_pp(p_ref_e_pp_indices) = d_perpendicular_e_pp_temp;

[d_perpendicular,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,p_test);
d_perpendicular = abs(d_perpendicular);

% PP Deviation from map CDF _______________________________________________

e_pp_cdf_error = 0:0.01:ceil(max(d_perpendicular)); % ensure to include all data by including maximum existing error
e_pp_cdf = [e_pp_cdf_error(2:end); zeros(1,length(e_pp_cdf_error)-1)];
e_pp_cdf(2,:) = histcounts(d_perpendicular,e_pp_cdf_error,'Normalization','cdf');

% Track-length deviation CDF ______________________________________________

time_cdf = sim_tgc.Time(data_indices(p_ref_pp_indices));
none_nan_data_selector = ~isnan(p_ref_pp(1,:));
d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp(:,none_nan_data_selector),1,2)].^2,1)));
p_test_temp = p_test(:,p_ref_pp_indices);
d_test = cumsum(sqrt(sum([zeros(2,1),diff(p_test_temp(:,none_nan_data_selector),1,2)].^2,1)));
% d_test = d_test(p_ref_pp_indices);
d_test = d_test - d_test(1);
e_track_length_orig = (d_p_ref_pp-d_test);
e_track_length = abs(e_track_length_orig);

e_tl_cdf_error = floor(min(e_track_length)):0.01:ceil(max(e_track_length)); % ensure to include all data by including maximum existing error
e_tl_cdf = [e_tl_cdf_error(2:end); zeros(1,length(e_tl_cdf_error)-1)];
e_tl_cdf(2,:) = histcounts(e_track_length,e_tl_cdf_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

% sim_tgc_error_data.e_time = time_orig;
sim_tgc_error_data.cdf_time = time_cdf;
sim_tgc_error_data.e_pp = d_perpendicular_e_pp;
sim_tgc_error_data.e_pp_cdf.MaxError_m = e_pp_cdf(1,:)';
sim_tgc_error_data.e_pp_cdf.Availability = e_pp_cdf(2,:)';
sim_tgc_error_data.e_tl = e_track_length_orig;
sim_tgc_error_data.e_tl_cdf.MaxError_m = e_tl_cdf(1,:)';
sim_tgc_error_data.e_tl_cdf.Availability = e_tl_cdf(2,:)';

fprintf('done!\n')

%% Optimized map

if exist('run_optimization_executed','var') && run_optimization_executed

    fprintf('\tCalculated map data...')

    % Perpendicular deviation from ref-map ________________________________

    %[~,p_test,~] = calcPerpendicularErrorToRefPoints(optim_map_abs_pos,p_test_optim_map_correspondence);
    p_test = optim_map_abs_pos;
    [d_perpendicular,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,p_test);
    d_perpendicular = abs(d_perpendicular);
    % d_perpendicular = d_perpendicular(~isnan(d_perpendicular));
    % p_ref_pp = p_ref_pp(:,~isnan(p_ref_pp(1,:)));

%     % Test ________________________________________________________________
%     p_test2 = optim_map_abs_pos;
%     [d_perpendicular2,p_ref_pp2,p_ref_pp_indices2] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,p_test2);
%     d_p_ref_pp2 = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp2,1,2)].^2,1)));
%     p_test_temp2 = p_test2(:,p_ref_pp_indices2);
%     d_test2 = cumsum(sqrt(sum([zeros(2,1),diff(p_test_temp2,1,2)].^2,1)));
%     e_track_length2 = (d_p_ref_pp2-d_test2);
%     e_track_length2 = abs(e_track_length2);
%     
%     plot(ref_map_abs_position_utm(1,:),ref_map_abs_position_utm(2,:)); hold all
%     plot(p_ref_pp2(1,:),p_ref_pp2(2,:));
%     plot(p_test2(1,:),p_test2(2,:));
%     axis equal   
    
    
    % PP Deviation from map CDF ___________________________________________

    e_pp_cdf_error = 0:0.01:ceil(max(d_perpendicular)); % ensure to include all data by including maximum existing error
    e_pp_cdf = [e_pp_cdf_error(2:end); zeros(1,length(e_pp_cdf_error)-1)];
    e_pp_cdf(2,:) = histcounts(d_perpendicular,e_pp_cdf_error,'Normalization','cdf');

    % Track-length deviation CDF __________________________________________

    d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp,1,2)].^2,1)));
    p_test_temp = p_test(:,p_ref_pp_indices);
    d_test = cumsum(sqrt(sum([zeros(2,1),diff(p_test_temp,1,2)].^2,1)));
    e_track_length_orig = (d_p_ref_pp-d_test);
    e_track_length = abs(e_track_length_orig);

    e_tl_cdf_error = floor(min(e_track_length)):0.01:ceil(max(e_track_length)); % ensure to include all data by including maximum existing error
    e_tl_cdf = [e_tl_cdf_error(2:end); zeros(1,length(e_tl_cdf_error)-1)];
    e_tl_cdf(2,:) = histcounts(e_track_length,e_tl_cdf_error,'Normalization','cdf');

    % Write data to struct_________________________________________________

    optim_map_error_data.e_pp = d_perpendicular;
    optim_map_error_data.e_pp_cdf.MaxError_m = e_pp_cdf(1,:)';
    optim_map_error_data.e_pp_cdf.Availability = e_pp_cdf(2,:)';
    optim_map_error_data.e_tl = e_track_length_orig;
    optim_map_error_data.e_tl_cdf.MaxError_m = e_tl_cdf(1,:)';
    optim_map_error_data.e_tl_cdf.Availability = e_tl_cdf(2,:)';

    fprintf('done!\n')
else
    
    optim_map_error_data.e_pp = nan;
    optim_map_error_data.e_pp_cdf.MaxError_m = nan;
    optim_map_error_data.e_pp_cdf.Availability = nan;
    optim_map_error_data.e_tl = nan;
    optim_map_error_data.e_tl_cdf.MaxError_m = nan;
    optim_map_error_data.e_tl_cdf.Availability = nan;
    
end % if

%% Save

% Update simout data ______________________________________________________

prompt_str = 'Do you want to save the updated processed simout data (j/n)?';
user_str = input(prompt_str,'s');

switch user_str
    case {'j','y','yes','ja'}
        simout_data_saved = false;
        save_sim_out_data = true;
    case {'n','no','nein'}
        save_sim_out_data = false;
    otherwise
        save_sim_out_data = false;
end % switch

if save_sim_out_data
    saveSimOutData
end % if

if exist('run_optimization_executed','var') && run_optimization_executed
    
    % Update optimout data ________________________________________________

    prompt_str = 'Do you want to update the processed optimout data (j/n)?';
    user_str = input(prompt_str,'s');

    switch user_str
        case {'j','y','yes','ja'}
            optimout_data_saved = false;
            save_optim_out_data = true;
        case {'n','no','nein'}
            save_optim_out_data = false;
        otherwise
            save_optim_out_data = false;
    end % switch

    if save_optim_out_data
        saveOptimData
    end %if

end % if

%% Helper Functions

function TF = isInArea(p_test,x_min,x_max,y_min,y_max)

if size(p_test,1) ~= 2
    error('calcDeviationToRefMap: Wrong dimension!')
end % if

TF = (p_test(1,:) > x_min) & ...
     (p_test(1,:) < x_max) & ...
     (p_test(2,:) > y_min) & ...
     (p_test(2,:) < y_max);

end % function
