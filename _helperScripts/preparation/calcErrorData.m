% 
% Calculate plot data for error-ellipses
%   
%   Author: Hanno Winter
%   Date: 13-Apr-2021; Last revision: 19-Jul-2021

%% Settings

% see: prepareSimOutData

%% Init

% if ~exist('error_ellipse_data_prepared','var') || error_ellipse_data_prepared == false
%     error_ellipse_data_prepared = false;
% elseif error_ellipse_data_prepared
%     fprintf('Already prepared error-ellipse data\n');
%     return  
% end % if

z_sigma_gain = sqrt(chi2inv(error_ellipse_confidence,2));

v_cdf_min = v_min_tgc_speed-v_min_tgc_hysteresis;
v_min_selector = abs(sim_imm.VelocityVehicle_ms) > v_cdf_min;
v_min_gps_selector = v_min_selector(simout_gnss_valid_flag.Data(:,1)==1);

% v_min_gps_selector = logical(simout_gnss_valid_flag.Data(:,1)) & v_min_selector;

%% GNSS

% Error-Ellipse ___________________________________________________________

gnss_error_ellipses_indices = 1:ceil(ts_error_ellipses/mean(diff(sim_gnss_filter_input.Time))):length(sim_gnss_filter_input.Time);
gnss_error_ellipses_indices = gnss_error_ellipses_indices(sim_gnss_filter_input.ErrorEllipseMajor_m(gnss_error_ellipses_indices)' < max_major_error_ellipses);
[~,gnss_error_ellipse] = ... 
    getLatLonErrorEllipsePoints( ... 
                                 sim_gnss_filter_input{gnss_error_ellipses_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                 sim_gnss_filter_input.ErrorEllipseMajor_m(gnss_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_gnss_filter_input.ErrorEllipseMinor_m(gnss_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_gnss_filter_input.ErrorEllipseOrientation_deg(gnss_error_ellipses_indices), ... 
                                 360 ... 
                               );

% Availability ____________________________________________________________

expected_gnss_time = simulink_time(1):natural_gnss_sample_time:simulink_time(end);

num_expected_gnss = simulink_time(end)/natural_gnss_sample_time+1; % +1 to include first time step at t=0
num_valid_gnss = length(sim_gnss_filter_input.Time);
num_missing_gnss_z = num_expected_gnss-num_valid_gnss;
sim_gnss_availability = num_valid_gnss/num_expected_gnss; 

% CDF _____________________________________________________________________

theta_along_track = simplifyHeadingD(90-sim_gnss_filter_input.Heading_deg,'180');
theta_cross_track = simplifyHeadingD(theta_along_track+90,'180');
z_sigma_at = nan(size(theta_along_track));
z_sigma_ct = nan(size(theta_along_track));
for i = 1:length(sim_gnss_filter_input.PositionCov)
    [z_sigma_at(i,1),~,~] = calcStandardDeviationCurve(sim_gnss_filter_input.PositionCov{i},theta_along_track(i),z_sigma_gain*1);
    [z_sigma_ct(i,1),~,~] = calcStandardDeviationCurve(sim_gnss_filter_input.PositionCov{i},theta_cross_track(i),z_sigma_gain*1);
end % for i

if exclude_stillstand_cdf
    num_missing_gnss_z = ceil(sum(v_min_gps_selector)*(1-sim_gnss_availability));
    z_sigma_major = sim_gnss_filter_input.ErrorEllipseMajor_m(v_min_gps_selector)*z_sigma_gain;
    z_sigma_major = [z_sigma_major;nan(num_missing_gnss_z,1)]; % add nan data to represent availability issues
    z_sigma_at_stripped = [z_sigma_at(v_min_gps_selector);nan(num_missing_gnss_z,1)]; % add nan data to represent availability issues
    z_sigma_ct_stripped = [z_sigma_ct(v_min_gps_selector);nan(num_missing_gnss_z,1)]; % add nan data to represent availability issues
else
    z_sigma_major = [z_sigma_major;nan(num_missing_gnss_z,1)]; % add nan data to represent availability issues
    z_sigma_at_stripped = [z_sigma_at;nan(num_missing_gnss_z,1)]; % add nan data to represent availability issues
    z_sigma_ct_stripped = [z_sigma_ct;nan(num_missing_gnss_z,1)]; % add nan data to represent availability issues
end % if

cdf_max_error = 0:0.01:ceil(max(z_sigma_major)); % ensure to include all data by including maximum existing error
sim_gnss_max_cdf = [cdf_max_error(2:end); zeros(1,length(cdf_max_error)-1)];
sim_gnss_max_cdf(2,:) = histcounts(z_sigma_major,cdf_max_error,'Normalization','cdf');

cdf_at_error = 0:0.01:ceil(max(z_sigma_at_stripped)); % ensure to include all data by including maximum existing error
sim_gnss_at_cdf = [cdf_at_error(2:end); zeros(1,length(cdf_at_error)-1)];
sim_gnss_at_cdf(2,:) = histcounts(z_sigma_at_stripped,cdf_at_error,'Normalization','cdf');

cdf_ct_error = 0:0.01:ceil(max(z_sigma_ct_stripped)); % ensure to include all data by including maximum existing error
sim_gnss_ct_cdf = [cdf_ct_error(2:end); zeros(1,length(cdf_ct_error)-1)];
sim_gnss_ct_cdf(2,:) = histcounts(z_sigma_ct_stripped,cdf_ct_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

sim_gnss_filter_input_error_data.e_time = sim_gnss_filter_input.Time;
sim_gnss_filter_input_error_data.s.sigma_at = z_sigma_at;
sim_gnss_filter_input_error_data.s.sigma_ct = z_sigma_ct;

sim_gnss_filter_input_error_data.ee.Latitude_deg = gnss_error_ellipse(1,:)';
sim_gnss_filter_input_error_data.ee.Longitude_deg = gnss_error_ellipse(2,:)';
sim_gnss_filter_input_error_data.ee.confidence = error_ellipse_confidence;
sim_gnss_filter_input_error_data.ee.sigma_gain = z_sigma_gain;

sim_gnss_filter_input_error_data.availability = sim_gnss_availability;

sim_gnss_filter_input_error_data.cdf_time = sim_gnss_filter_input.Time(v_min_gps_selector);

sim_gnss_filter_input_error_data.cdf_max.MaxError_m = sim_gnss_max_cdf(1,:)';
sim_gnss_filter_input_error_data.cdf_max.Availability = sim_gnss_max_cdf(2,:)';
sim_gnss_filter_input_error_data.cdf_max.confidence = error_ellipse_confidence;
sim_gnss_filter_input_error_data.cdf_max.sigma_gain = z_sigma_gain;

sim_gnss_filter_input_error_data.cdf_at.MaxError_m = sim_gnss_at_cdf(1,:)';
sim_gnss_filter_input_error_data.cdf_at.Availability = sim_gnss_at_cdf(2,:)';
sim_gnss_filter_input_error_data.cdf_at.confidence = error_ellipse_confidence;
sim_gnss_filter_input_error_data.cdf_at.sigma_gain = z_sigma_gain;

sim_gnss_filter_input_error_data.cdf_ct.MaxError_m = sim_gnss_ct_cdf(1,:)';
sim_gnss_filter_input_error_data.cdf_ct.Availability = sim_gnss_ct_cdf(2,:)';
sim_gnss_filter_input_error_data.cdf_ct.confidence = error_ellipse_confidence;
sim_gnss_filter_input_error_data.cdf_ct.sigma_gain = z_sigma_gain;

%% EKF

% Error-Ellipse ___________________________________________________________

ekf_error_ellipses_indices = 1:ceil(ts_error_ellipses/mean(diff(sim_ekf.Time))):length(sim_ekf.Time);
ekf_error_ellipses_indices = ekf_error_ellipses_indices(sim_ekf.ErrorEllipseMajor_m(ekf_error_ellipses_indices)' < max_major_error_ellipses);
[~,ekf_error_ellipse] = ... 
    getLatLonErrorEllipsePoints( ... 
                                 sim_ekf{ekf_error_ellipses_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                 sim_ekf.ErrorEllipseMajor_m(ekf_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_ekf.ErrorEllipseMinor_m(ekf_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_ekf.ErrorEllipseOrientation_deg(ekf_error_ellipses_indices), ... 
                                 360 ... 
                               );

% CDF _____________________________________________________________________

theta_along_track = simplifyHeadingD(90-sim_ekf.Heading_deg,'180');
theta_cross_track = simplifyHeadingD(theta_along_track+90,'180');
z_sigma_at = nan(size(theta_along_track));
z_sigma_ct = nan(size(theta_along_track));
for i = 1:length(sim_ekf.PositionCov)
    [z_sigma_at(i,1),~,~] = calcStandardDeviationCurve(sim_ekf.PositionCov{i},theta_along_track(i),z_sigma_gain*1);
    [z_sigma_ct(i,1),~,~] = calcStandardDeviationCurve(sim_ekf.PositionCov{i},theta_cross_track(i),z_sigma_gain*1);
end % for i
z_sigma_major = sim_ekf.ErrorEllipseMajor_m*z_sigma_gain;

expected_gnss_time = simulink_time(1):natural_gnss_sample_time:simulink_time(end);
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
z_sigma_major = z_sigma_major(data_selector);
z_sigma_at_stripped = z_sigma_at(data_selector);
z_sigma_ct_stripped = z_sigma_ct(data_selector);

cdf_max_error = 0:0.01:ceil(max(z_sigma_major)); % ensure to include all data by including maximum existing error
sim_ekf_max_cdf = [cdf_max_error(2:end); zeros(1,length(cdf_max_error)-1)];
sim_ekf_max_cdf(2,:) = histcounts(z_sigma_major,cdf_max_error,'Normalization','cdf');

cdf_at_error = 0:0.01:ceil(max(z_sigma_at_stripped)); % ensure to include all data by including maximum existing error
sim_ekf_at_cdf = [cdf_at_error(2:end); zeros(1,length(cdf_at_error)-1)];
sim_ekf_at_cdf(2,:) = histcounts(z_sigma_at_stripped,cdf_at_error,'Normalization','cdf');

cdf_ct_error = 0:0.01:ceil(max(z_sigma_ct_stripped)); % ensure to include all data by including maximum existing error
sim_ekf_ct_cdf = [cdf_ct_error(2:end); zeros(1,length(cdf_ct_error)-1)];
sim_ekf_ct_cdf(2,:) = histcounts(z_sigma_ct_stripped,cdf_ct_error,'Normalization','cdf');

% Write data to struct ____________________________________________________

sim_ekf_error_data.e_time = sim_ekf.Time;
sim_ekf_error_data.s.sigma_at = z_sigma_at;
sim_ekf_error_data.s.sigma_ct = z_sigma_ct;

sim_ekf_error_data.ee.Latitude_deg = ekf_error_ellipse(1,:)';
sim_ekf_error_data.ee.Longitude_deg = ekf_error_ellipse(2,:)';
sim_ekf_error_data.ee.confidence = error_ellipse_confidence;
sim_ekf_error_data.ee.sigma_gain = z_sigma_gain;

sim_ekf_error_data.cdf_time = sim_ekf.Time(data_selector);

sim_ekf_error_data.cdf_max.MaxError_m = sim_ekf_max_cdf(1,:)';
sim_ekf_error_data.cdf_max.Availability = sim_ekf_max_cdf(2,:)';
sim_ekf_error_data.cdf_max.confidence = error_ellipse_confidence;
sim_ekf_error_data.cdf_max.sigma_gain = z_sigma_gain;

sim_ekf_error_data.cdf_at.MaxError_m = sim_ekf_at_cdf(1,:)';
sim_ekf_error_data.cdf_at.Availability = sim_ekf_at_cdf(2,:)';
sim_ekf_error_data.cdf_at.confidence = error_ellipse_confidence;
sim_ekf_error_data.cdf_at.sigma_gain = z_sigma_gain;

sim_ekf_error_data.cdf_ct.MaxError_m = sim_ekf_ct_cdf(1,:)';
sim_ekf_error_data.cdf_ct.Availability = sim_ekf_ct_cdf(2,:)';
sim_ekf_error_data.cdf_ct.confidence = error_ellipse_confidence;
sim_ekf_error_data.cdf_ct.sigma_gain = z_sigma_gain;

%% IMM

% Error-Ellipse ___________________________________________________________

imm_error_ellipses_indices = 1:ceil(ts_error_ellipses/mean(diff(sim_imm.Time))):length(sim_imm.Time);
imm_error_ellipses_indices = imm_error_ellipses_indices(sim_imm.ErrorEllipseMajor_m(imm_error_ellipses_indices)' < max_major_error_ellipses);
[~,imm_error_ellipse] = ... 
    getLatLonErrorEllipsePoints( ... 
                                 sim_imm{imm_error_ellipses_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                 sim_imm.ErrorEllipseMajor_m(imm_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_imm.ErrorEllipseMinor_m(imm_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_imm.ErrorEllipseOrientation_deg(imm_error_ellipses_indices), ... 
                                 360 ... 
                               );

% CDF _____________________________________________________________________

% data_selector = ismember(sim_imm.Time,expected_gnss_time);
% data_selector = true(size(sim_imm.Time));

theta_along_track = simplifyHeadingD(90-sim_imm.Heading_deg,'180');
theta_cross_track = simplifyHeadingD(theta_along_track+90,'180');
z_sigma_at = nan(size(theta_along_track));
z_sigma_ct = nan(size(theta_along_track));
for i = 1:length(sim_imm.PositionCov)
    [z_sigma_at(i,1),~,~] = calcStandardDeviationCurve(sim_imm.PositionCov{i},theta_along_track(i),z_sigma_gain*1);
    [z_sigma_ct(i,1),~,~] = calcStandardDeviationCurve(sim_imm.PositionCov{i},theta_cross_track(i),z_sigma_gain*1);
end % for i
z_sigma_major = sim_imm.ErrorEllipseMajor_m*z_sigma_gain;

time_selector = ismember(sim_imm.Time,expected_gnss_time);
if exclude_stillstand_cdf && harmonize_with_gps
    data_selector = time_selector & v_min_selector;
elseif exclude_stillstand_cdf
    data_selector = v_min_selector;
elseif harmonize_with_gps
    data_selector = time_selector;
else
    data_selector = true(size(v_min_selector));
end % if
z_sigma_major = z_sigma_major(data_selector);
z_sigma_at_stripped = z_sigma_at(data_selector);
z_sigma_ct_stripped = z_sigma_ct(data_selector);

cdf_max_error = 0:0.01:ceil(max(z_sigma_major)); % ensure to include all data by including maximum existing error
sim_imm_max_cdf = [cdf_max_error(2:end); zeros(1,length(cdf_max_error)-1)];
sim_imm_max_cdf(2,:) = histcounts(z_sigma_major,cdf_max_error,'Normalization','cdf');

cdf_at_error = 0:0.01:ceil(max(z_sigma_at_stripped)); % ensure to include all data by including maximum existing error
sim_imm_at_cdf = [cdf_at_error(2:end); zeros(1,length(cdf_at_error)-1)];
sim_imm_at_cdf(2,:) = histcounts(z_sigma_at_stripped,cdf_at_error,'Normalization','cdf');

cdf_ct_error = 0:0.01:ceil(max(z_sigma_ct_stripped)); % ensure to include all data by including maximum existing error
sim_imm_ct_cdf = [cdf_ct_error(2:end); zeros(1,length(cdf_ct_error)-1)];
sim_imm_ct_cdf(2,:) = histcounts(z_sigma_ct_stripped,cdf_ct_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

sim_imm_error_data.e_time = sim_imm.Time;
sim_imm_error_data.s.sigma_at = z_sigma_at;
sim_imm_error_data.s.sigma_ct = z_sigma_ct;

sim_imm_error_data.ee.Latitude_deg = imm_error_ellipse(1,:)';
sim_imm_error_data.ee.Longitude_deg = imm_error_ellipse(2,:)';
sim_imm_error_data.ee.confidence = error_ellipse_confidence;
sim_imm_error_data.ee.sigma_gain = z_sigma_gain;

sim_imm_error_data.cdf_time = sim_imm.Time(data_selector);

sim_imm_error_data.cdf_max.MaxError_m = sim_imm_max_cdf(1,:)';
sim_imm_error_data.cdf_max.Availability = sim_imm_max_cdf(2,:)';
sim_imm_error_data.cdf_max.confidence = error_ellipse_confidence;
sim_imm_error_data.cdf_max.sigma_gain = z_sigma_gain;

sim_imm_error_data.cdf_at.MaxError_m = sim_imm_at_cdf(1,:)';
sim_imm_error_data.cdf_at.Availability = sim_imm_at_cdf(2,:)';
sim_imm_error_data.cdf_at.confidence = error_ellipse_confidence;
sim_imm_error_data.cdf_at.sigma_gain = z_sigma_gain;

sim_imm_error_data.cdf_ct.MaxError_m = sim_imm_ct_cdf(1,:)';
sim_imm_error_data.cdf_ct.Availability = sim_imm_ct_cdf(2,:)';
sim_imm_error_data.cdf_ct.confidence = error_ellipse_confidence;
sim_imm_error_data.cdf_ct.sigma_gain = z_sigma_gain;

%% TGC

% Error-Ellipse ___________________________________________________________

tgc_error_ellipses_indices = 1:ceil(ts_error_ellipses/mean(diff(sim_tgc.Time))):length(sim_tgc.Time);
tgc_error_ellipses_indices = tgc_error_ellipses_indices(sim_tgc.ErrorEllipseMajor_m(tgc_error_ellipses_indices)' < max_major_error_ellipses);
[~,tgc_error_ellipse] = ... 
    getLatLonErrorEllipsePoints( ... 
                                 sim_tgc{tgc_error_ellipses_indices,{'Latitude_deg','Longitude_deg'}}', ... 
                                 sim_tgc.ErrorEllipseMajor_m(tgc_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_tgc.ErrorEllipseMinor_m(tgc_error_ellipses_indices)*z_sigma_gain, ... 
                                 sim_tgc.ErrorEllipseOrientation_deg(tgc_error_ellipses_indices), ... 
                                 360 ... 
                               );
% CDF _____________________________________________________________________

theta_along_track = simplifyHeadingD(90-sim_tgc.Heading_deg,'180');
theta_cross_track = simplifyHeadingD(theta_along_track+90,'180');
z_sigma_at = nan(size(theta_along_track));
z_sigma_ct = nan(size(theta_along_track));
for i = 1:length(sim_tgc.PositionCov)
    [z_sigma_at(i,1),~,~] = calcStandardDeviationCurve(sim_tgc.PositionCov{i},theta_along_track(i),z_sigma_gain*1);
    [z_sigma_ct(i,1),~,~] = calcStandardDeviationCurve(sim_tgc.PositionCov{i},theta_cross_track(i),z_sigma_gain*1);
end % for i
z_sigma_major = sim_tgc.ErrorEllipseMajor_m*z_sigma_gain;

time_selector = ismember(sim_tgc.Time,expected_gnss_time);
if exclude_stillstand_cdf && harmonize_with_gps
    data_selector = time_selector & v_min_selector;
elseif exclude_stillstand_cdf
    data_selector = v_min_selector;
elseif harmonize_with_gps
    data_selector = time_selector;
else
    data_selector = true(size(v_min_selector));
end % if
z_sigma_major = z_sigma_major(data_selector);
z_sigma_at_stripped = z_sigma_at(data_selector);
z_sigma_ct_stripped = z_sigma_ct(data_selector);

cdf_max_error = 0:0.01:ceil(max(z_sigma_major)); % ensure to include all data by including maximum existing error
sim_tgc_max_cdf = [cdf_max_error(2:end); zeros(1,length(cdf_max_error)-1)];
sim_tgc_max_cdf(2,:) = histcounts(z_sigma_major,cdf_max_error,'Normalization','cdf');

cdf_at_error = 0:0.01:ceil(max(z_sigma_at_stripped)); % ensure to include all data by including maximum existing error
sim_tgc_at_cdf = [cdf_at_error(2:end); zeros(1,length(cdf_at_error)-1)];
sim_tgc_at_cdf(2,:) = histcounts(z_sigma_at_stripped,cdf_at_error,'Normalization','cdf');

cdf_ct_error = 0:0.01:ceil(max(z_sigma_ct_stripped)); % ensure to include all data by including maximum existing error
sim_tgc_ct_cdf = [cdf_ct_error(2:end); zeros(1,length(cdf_ct_error)-1)];
sim_tgc_ct_cdf(2,:) = histcounts(z_sigma_ct_stripped,cdf_ct_error,'Normalization','cdf');

% Write data to struct_____________________________________________________

sim_tgc_error_data.e_time = sim_tgc.Time;
sim_tgc_error_data.s.sigma_at = z_sigma_at;
sim_tgc_error_data.s.sigma_ct = z_sigma_ct;

sim_tgc_error_data.ee.Latitude_deg = tgc_error_ellipse(1,:)';
sim_tgc_error_data.ee.Longitude_deg = tgc_error_ellipse(2,:)';
sim_tgc_error_data.ee.confidence = error_ellipse_confidence;
sim_tgc_error_data.ee.sigma_gain = z_sigma_gain;

sim_tgc_error_data.cdf_time = sim_tgc.Time(data_selector);

sim_tgc_error_data.cdf_max.MaxError_m = sim_tgc_max_cdf(1,:)';
sim_tgc_error_data.cdf_max.Availability = sim_tgc_max_cdf(2,:)';
sim_tgc_error_data.cdf_max.confidence = error_ellipse_confidence;
sim_tgc_error_data.cdf_max.sigma_gain = z_sigma_gain;

sim_tgc_error_data.cdf_at.MaxError_m = sim_tgc_at_cdf(1,:)';
sim_tgc_error_data.cdf_at.Availability = sim_tgc_at_cdf(2,:)';
sim_tgc_error_data.cdf_at.confidence = error_ellipse_confidence;
sim_tgc_error_data.cdf_at.sigma_gain = z_sigma_gain;

sim_tgc_error_data.cdf_ct.MaxError_m = sim_tgc_ct_cdf(1,:)';
sim_tgc_error_data.cdf_ct.Availability = sim_tgc_ct_cdf(2,:)';
sim_tgc_error_data.cdf_ct.confidence = error_ellipse_confidence;
sim_tgc_error_data.cdf_ct.sigma_gain = z_sigma_gain;

%% Uncertainty Growth Rate

v_gps_out_min = v_min_tgc_speed-v_min_tgc_hysteresis;
gps_out_v_min_indices = find(abs(sim_imm.VelocityVehicle_ms) > v_gps_out_min);

% Find GNSS outages start and end indices _________________________________

gnss_availability_str = sprintf('%d',simout_gnss_valid_flag.Data(:,1));

gnss_out_start_indices = strfind(gnss_availability_str,'10')'+1;
gnss_out_end_indices = strfind(gnss_availability_str,'01')';
if strcmp(gnss_availability_str(1),'0')
    gnss_out_start_indices = [1;gnss_out_start_indices];
end % if
if strcmp(gnss_availability_str(end),'0')
    gnss_out_end_indices = [gnss_out_end_indices;length(gnss_availability_str)];
end % if

% only use indices where vehicle is moving
gps_out_v_min_selector = (ismember(gnss_out_start_indices,gps_out_v_min_indices) & ismember(gnss_out_end_indices,gps_out_v_min_indices));
gnss_out_start_indices = gnss_out_start_indices(gps_out_v_min_selector);
gnss_out_end_indices = gnss_out_end_indices(gps_out_v_min_selector);
gnss_out_indices = [gnss_out_start_indices, gnss_out_end_indices];

% Get outage times and duration ___________________________________________

gnss_out_start_times = simout_gnss_valid_flag.Time(gnss_out_start_indices);
gnss_out_end_times = simout_gnss_valid_flag.Time(gnss_out_end_indices);
gnss_out_times = [gnss_out_start_times,gnss_out_end_times];
gnss_out_diff_times = diff(gnss_out_times,1,2);

% Get uncertainties and uncertainty differences ___________________________

gnss_out_ekf_sigma_start = sim_ekf_error_data.s.sigma_ct(gnss_out_start_indices);
gnss_out_ekf_sigma_end = sim_ekf_error_data.s.sigma_ct(gnss_out_end_indices);
gnss_out_ekf_sigma = [gnss_out_ekf_sigma_start, gnss_out_ekf_sigma_end];
gnss_out_ekf_delta_sigma = diff(gnss_out_ekf_sigma,1,2);

gnss_out_imm_sigma_start = sim_imm_error_data.s.sigma_ct(gnss_out_start_indices);
gnss_out_imm_sigma_end = sim_imm_error_data.s.sigma_ct(gnss_out_end_indices);
gnss_out_imm_sigma = [gnss_out_imm_sigma_start, gnss_out_imm_sigma_end];
gnss_out_imm_delta_sigma = diff(gnss_out_imm_sigma,1,2);

gnss_out_tgc_sigma_start = sim_tgc_error_data.s.sigma_ct(gnss_out_start_indices);
gnss_out_tgc_sigma_end = sim_tgc_error_data.s.sigma_ct(gnss_out_end_indices);
gnss_out_tgc_sigma = [gnss_out_tgc_sigma_start, gnss_out_tgc_sigma_end];
gnss_out_tgc_delta_sigma = diff(gnss_out_tgc_sigma,1,2);

% Calculate uncertainty growth ____________________________________________

outage_selector = (gnss_out_diff_times > 1);

sim_ekf_error_data.sigma_growth = gnss_out_ekf_delta_sigma(outage_selector)./gnss_out_diff_times(outage_selector);
sim_ekf_error_data.sigma_growth = sim_ekf_error_data.sigma_growth(sim_ekf_error_data.sigma_growth>0);
sim_imm_error_data.sigma_growth = gnss_out_imm_delta_sigma(outage_selector)./gnss_out_diff_times(outage_selector);
sim_imm_error_data.sigma_growth = sim_imm_error_data.sigma_growth(sim_imm_error_data.sigma_growth>0);
sim_tgc_error_data.sigma_growth = gnss_out_tgc_delta_sigma(outage_selector)./gnss_out_diff_times(outage_selector);
sim_tgc_error_data.sigma_growth = sim_tgc_error_data.sigma_growth(sim_tgc_error_data.sigma_growth>0);

% close all
% figure;histogram(sim_ekf_error_data.sigma_growth,-10:0.1:10);
% figure;histogram(sim_imm_error_data.sigma_growth,-10:0.1:10);
% figure;histogram(sim_tgc_error_data.sigma_growth,-10:0.1:10);

% [mean(sim_ekf_error_data.sigma_growth), median(sim_ekf_error_data.sigma_growth)]
% [mean(sim_imm_error_data.sigma_growth), median(sim_imm_error_data.sigma_growth)]
% [mean(sim_tgc_error_data.sigma_growth), median(sim_tgc_error_data.sigma_growth)]

%% Finish

% error_ellipse_data_prepared = true;
