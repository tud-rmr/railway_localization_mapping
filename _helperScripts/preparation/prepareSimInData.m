% 
% Prepare data for simulink model
%   
%   Author: Hanno Winter
%   Date: 20-Mar-2021; Last revision: 18-Jul-2021

%% Init

if ~exist('init_localization_executed','var')
    warning('prepareSimInputData.m: This script is supposed to be called from ''initLocalization.m''!')
    return
else
    prepare_sim_input_data_executed = false;
end % if

%% Prepare data

% Shrink data to single timetable _________________________________________

imu_data_selector = cellfun(@(cell) size(cell,1)>1,imu_data);
if sum(imu_data_selector) > 1
    error('initLocalization.m: Session based data not supported!');
else
    imu_data = imu_data{imu_data_selector};
end % if

gnss_data_selector = cellfun(@(cell) size(cell,1)>1,gnss_data);
if sum(imu_data_selector) > 1
    error('initLocalization.m: Session based data not supported!');
else
    gnss_data = gnss_data{gnss_data_selector};
end % if

% Find first valid GNSS measurement _______________________________________

gnss_data_validity = ~isnan(gnss_data.Latitude_deg) & ...
                     ~isnan(gnss_data.Longitude_deg) & ...
                     ~isnan(gnss_data.AltitudeEllipsoid_m) & ...
                     ~isnan(gnss_data.VelocityNorth_ms) & ...
                     ~isnan(gnss_data.VelocityEast_ms) & ...
                     ~isnan(gnss_data.VelocityDown_ms) & ...
                     ~isnan(gnss_data.Heading_deg);
                 
if any(~isnan(gnss_data.PDOP)) && ~isempty(max_pdop)
    pdop_selector = ~isnan(gnss_data.PDOP);
    gnss_data_validity(pdop_selector,1) = gnss_data_validity(pdop_selector,1) & ...
                                          (gnss_data.PDOP(pdop_selector) < max_pdop);
end % if

first_valid_gnss_index = find(gnss_data_validity==true,1);
first_valid_imu_idnex = max(find(imu_data.Time>=gnss_data.Time(first_valid_gnss_index),1)-1,1);

% Trim input data to first valid GNSS signal ______________________________

imu_data = imu_data(first_valid_imu_idnex:end,:);
gnss_data = gnss_data(first_valid_gnss_index:end,:);
gnss_data_validity = gnss_data_validity(first_valid_gnss_index:end);

%% Set sample rates and time vectors 

time_imu = seconds(imu_data.Time);
time_gnss = seconds(gnss_data.Time);

t0 = max([time_imu(1),time_gnss(1)]);
time_imu = time_imu - t0;
time_gnss = time_gnss - t0;

sample_time = min([imu_sample_time,gnss_sample_time]);

[N,edges] = histcounts(diff(time_imu));
[~,max_index] = max(N);
natural_imu_sample_time = mean(edges([max_index,max_index+1]));
natural_imu_sample_time = round(natural_imu_sample_time*1e3)/1e3;

[N,edges] = histcounts(diff(time_gnss));
[~,max_index] = max(N);
natural_gnss_sample_time = mean(edges([max_index,max_index+1]));
natural_gnss_sample_time = round(natural_gnss_sample_time*1e1)/1e1;

natural_sample_time = min([natural_imu_sample_time,natural_gnss_sample_time]);

if natural_sample_time > sample_time
    sample_time = natural_sample_time;
    warning('initLocalization.m: Reset simulation sample time to natural rate!');
end % if

if natural_imu_sample_time > imu_sample_time
    imu_sample_time = natural_imu_sample_time;
    warning('initLocalization.m: Reset IMU sample time to natural rate!');
end % if

if natural_gnss_sample_time > gnss_sample_time
    gnss_sample_time = natural_gnss_sample_time;
    warning('initLocalization.m: Reset GNSS sample time to natural rate!');
end % if

imu_utc_time = datetime(seconds(imu_data.Time),'ConvertFrom','posixtime');
gnss_utc_time = datetime(seconds(gnss_data.Time),'ConvertFrom','posixtime');

%% Insert artificial GNSS outages

gnss_data_validity_original = gnss_data_validity; % Save original validity data
insertArtificialOutages

%% Set variables for simulink

imu_lever_arm = imu_parameters{1,{'LeverArmX_m','LeverArmY_m','LeverArmZ_m'}};
gnss_lever_arm = gnss_parameters{1,{'LeverArmX_m','LeverArmY_m','LeverArmZ_m'}};

acc_data = imu_data{:,{'AccX_mss','AccY_mss','AccZ_mss'}};
acc_data = timeseries(fillmissing(acc_data,'previous'),time_imu');
acc_data_validity = timeseries(1,0);

gyro_data = imu_data{:,{'TurnRateX_degs','TurnRateY_degs','TurnRateZ_degs'}};
gyro_data = timeseries(fillmissing(gyro_data,'previous'),time_imu');
gyro_data_validity = timeseries(1,0);

gnss_sim_data = gnss_data{:,{'Latitude_deg','Longitude_deg','AltitudeEllipsoid_m','VelocityNorth_ms','VelocityEast_ms','VelocityDown_ms','Heading_deg'}};
gnss_sim_data = timeseries(fillmissing(gnss_sim_data,'previous'),time_gnss');
gnss_sim_data_validity = timeseries(gnss_data_validity,time_gnss');

gnss_uncertainty_sim_data = gnss_data{:,{'LatitudeSigma_m','LongitudeSigma_m','VelocityGroundSigma_ms','HeadingSigma_deg'}};
gnss_uncertainty_sim_data = timeseries(fillmissing(gnss_uncertainty_sim_data,'previous'),time_gnss');

v_min_tgc = [v_min_tgc_speed;v_min_tgc_hysteresis];

%% Finish script

prepare_sim_input_data_executed = true;
