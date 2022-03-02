% 
% Prepare output data from 'localization.slx' for further usage
%   
%   Author: Hanno Winter
%   Date: 26-Mar-2021; Last revision: 14-Dec-2021

%% Settings

sim_map_density = 1; % Railway-map resolution in m

ts_error_ellipses = 1; % time between error-ellipse plots in sec
max_major_error_ellipses = 50; % maximum major semi-axis size of error ellipses to plot
error_ellipse_confidence = 0.99; % error ellipse confidence
exclude_stillstand_cdf = 1; % exclude stanstills in CDF calculation
harmonize_with_gps = 1; % harmonize data used for CDF to GPS sample size

%% Init

if ~exist('simout_data_prepared','var') || simout_data_prepared == false
    simout_data_prepared = false;
    fprintf('Preparing Simulink output data:\n');
elseif simout_data_prepared
    fprintf('Already prepared Simulink output data\n');
    return  
end % if

simulink_time = simout_sensor_data.Time;
p_0_utm = simout_p_0_utm.Data(:);

num_processing_steps = 9;

%% Strapdown Data

fprintf('\t 1/%i ''Strapdown data''...',num_processing_steps);

% IMU _____________________________________________________________________
sim_a_ib_b = table();
sim_a_ib_b.Time = simout_a_ib_b.Time;
sim_a_ib_b.UtcTime = datetime(t0+sim_a_ib_b.Time,'ConvertFrom','posixtime');
sim_a_ib_b.AccX_mss = simout_a_ib_b.Data(:,1);
sim_a_ib_b.AccY_mss = simout_a_ib_b.Data(:,2);
sim_a_ib_b.AccZ_mss = simout_a_ib_b.Data(:,3);

sim_w_ib_b = table();
sim_w_ib_b.Time = simout_w_ib_b.Time;
sim_w_ib_b.UtcTime = datetime(t0+sim_w_ib_b.Time,'ConvertFrom','posixtime');
sim_w_ib_b.TurnRateX_degs = simout_w_ib_b.Data(:,1);
sim_w_ib_b.TurnRateY_degs = simout_w_ib_b.Data(:,2);
sim_w_ib_b.TurnRateZ_degs = simout_w_ib_b.Data(:,3);

% GNSS ____________________________________________________________________
sim_gnss = table();
sim_gnss.Time = simout_gnss.gnss_sensor_data_corrected.Time;
sim_gnss.UtcTime = datetime(t0+sim_gnss.Time,'ConvertFrom','posixtime');
sim_gnss.Latitude_deg = simout_gnss.gnss_sensor_data_corrected.Data(:,1);
sim_gnss.Longitude_deg = simout_gnss.gnss_sensor_data_corrected.Data(:,2);
sim_gnss.AltitudeEllipsoid_m = simout_gnss.gnss_sensor_data_corrected.Data(:,3);
sim_gnss.VelocityNorth_ms = simout_gnss.gnss_sensor_data_corrected.Data(:,4);
sim_gnss.VelocityEast_ms = simout_gnss.gnss_sensor_data_corrected.Data(:,5);
sim_gnss.VelocityDown_ms = simout_gnss.gnss_sensor_data_corrected.Data(:,6);
sim_gnss.Heading_deg = simout_gnss.gnss_sensor_data_corrected.Data(:,7);

% INS _____________________________________________________________________

sim_attitude = table();
sim_attitude.Time = simout_attitude_deg.Time;
sim_attitude.UtcTime = datetime(t0+sim_attitude.Time,'ConvertFrom','posixtime');
sim_attitude.Roll_deg = simout_attitude_deg.Data(:,1);
sim_attitude.Pitch_deg = simout_attitude_deg.Data(:,2);
sim_attitude.Heading_deg = simout_attitude_deg.Data(:,3);

sim_attitude_raw = table();
sim_attitude_raw.Time = simout_attitude_raw_deg.Time;
sim_attitude_raw.UtcTime = datetime(t0+sim_attitude_raw.Time,'ConvertFrom','posixtime');
sim_attitude_raw.Roll_deg = simout_attitude_raw_deg.Data(:,1);
sim_attitude_raw.Pitch_deg = simout_attitude_raw_deg.Data(:,2);
sim_attitude_raw.Heading_deg = simout_attitude_raw_deg.Data(:,3);

sim_v_ned = table();
sim_v_ned.Time = simout_v_ned_n.Time;
sim_v_ned.UtcTime = datetime(t0+sim_v_ned.Time,'ConvertFrom','posixtime');
sim_v_ned.VelocityNorth_ms = simout_v_ned_n.Data(:,1);
sim_v_ned.VelocityEast_ms = simout_v_ned_n.Data(:,2);
sim_v_ned.VelocityDown_ms = simout_v_ned_n.Data(:,3);

sim_w_nb_b = table();
sim_w_nb_b.Time = simout_w_nb_b.Time;
sim_w_nb_b.UtcTime = datetime(t0+sim_w_nb_b.Time,'ConvertFrom','posixtime');
sim_w_nb_b.TurnRateX_degs = simout_w_nb_b.Data(:,1) * 180/pi;
sim_w_nb_b.TurnRateY_degs = simout_w_nb_b.Data(:,2) * 180/pi;
sim_w_nb_b.TurnRateZ_degs = simout_w_nb_b.Data(:,3) * 180/pi;

sim_a_eb_n = table();
sim_a_eb_n.Time = simout_a_eb_n.Time;
sim_a_eb_n.UtcTime = datetime(t0+sim_a_eb_n.Time,'ConvertFrom','posixtime');
sim_a_eb_n.AccX_mss = simout_a_eb_n.Data(:,1);
sim_a_eb_n.AccY_mss = simout_a_eb_n.Data(:,2);
sim_a_eb_n.AccZ_mss = simout_a_eb_n.Data(:,3);

sim_q_b_n = table();
sim_q_b_n.Time = simout_q_b_n.Time;
sim_q_b_n.UtcTime = datetime(t0+sim_q_b_n.Time,'ConvertFrom','posixtime');
sim_q_b_n.q_b_n = mat2cell(simout_q_b_n.Data',4,ones(1,size(simout_q_b_n.Data,1)))';

% Misc ____________________________________________________________________

sim_a_bias = table();
sim_a_bias.Time = simout_a_bias_b.Time;
sim_a_bias.UtcTime = datetime(t0+sim_a_bias.Time,'ConvertFrom','posixtime');
sim_a_bias.BiasAccX_mss = simout_a_bias_b.Data(:,1);
sim_a_bias.BiasAccY_mss = simout_a_bias_b.Data(:,2);
sim_a_bias.BiasAccZ_mss = simout_a_bias_b.Data(:,3);

sim_w_bias = table();
sim_w_bias.Time = simout_w_bias_b.Time;
sim_w_bias.UtcTime = datetime(t0+sim_w_bias.Time,'ConvertFrom','posixtime');
sim_w_bias.BiasTurnRateX_mss = simout_w_bias_b.Data(:,1);
sim_w_bias.BiasTurnRateY_mss = simout_w_bias_b.Data(:,2);
sim_w_bias.BiasTurnRateZ_mss = simout_w_bias_b.Data(:,3);

sim_standstill_flag = table();
sim_standstill_flag.Time = simout_stillstand_flag.Time;
sim_standstill_flag.UtcTime = datetime(t0+sim_standstill_flag.Time,'ConvertFrom','posixtime');
sim_standstill_flag.standstill_flag = squeeze(simout_stillstand_flag.Data);

sim_w_energy = table();
sim_w_energy.Time = simout_w_energy.Time;
sim_w_energy.UtcTime = datetime(t0+sim_w_energy.Time,'ConvertFrom','posixtime');
sim_w_energy.w_energy = simout_w_energy.Data;

sim_a_e = table();
sim_a_e.Time = simout_a_e.Time;
sim_a_e.UtcTime = datetime(t0+sim_a_e.Time,'ConvertFrom','posixtime');
sim_a_e.a_e = simout_a_e.Data;

fprintf('done!\n');

%% Filter Input Data

fprintf('\t 2/%i ''Filter input data''...',num_processing_steps);

% Calculations ____________________________________________________________
x_utm = p_0_utm(1)+simout_sensor_data.Data(:,1);
y_utm = p_0_utm(2)+simout_sensor_data.Data(:,2);
[lat,lon] = utm2ll(x_utm,y_utm,p_0_utm(3),'wgs84');

major_minor_confidence = chi2cdf(1^2,2); % calculate probability for 1x sigma 2D-error-ellipse
[std_major,std_minor,alpha] = calcErrorEllipseParameters(simout_sensor_data_P.Data(1:2,1:2,:),major_minor_confidence);

P_pos = reshape(simout_sensor_data_P.Data(1:2,1:2,:),2,2*length(simout_sensor_data_P.Time));
P_pos = mat2cell(P_pos,2,repmat(2,1,length(simout_sensor_data_P.Time)))';

gnss_valid_selector = logical(simout_gnss_valid_flag.Data(:,1));

% Write to table __________________________________________________________
sim_gnss_filter_input = table();
sim_gnss_filter_input.Time = simulink_time(gnss_valid_selector);
sim_gnss_filter_input.UtcTime = datetime(t0+sim_gnss_filter_input.Time,'ConvertFrom','posixtime');
sim_gnss_filter_input.UtmEast_m = x_utm(gnss_valid_selector);
sim_gnss_filter_input.UtmNorth_m = y_utm(gnss_valid_selector);
sim_gnss_filter_input.Latitude_deg = lat(gnss_valid_selector);
sim_gnss_filter_input.Longitude_deg = lon(gnss_valid_selector);
sim_gnss_filter_input.VelocityGround_ms = simout_sensor_data.Data(gnss_valid_selector,3);
sim_gnss_filter_input.Heading_deg = simplifyHeadingD(rad2deg(simout_sensor_data.Data(gnss_valid_selector,4)),'180');
sim_gnss_filter_input.LatitudeSigma_m = squeeze(simout_sensor_data_P.Data(2,2,gnss_valid_selector));
sim_gnss_filter_input.LongitudeSigma_m = squeeze(simout_sensor_data_P.Data(1,1,gnss_valid_selector));
sim_gnss_filter_input.PositionCov = P_pos(gnss_valid_selector,:);
sim_gnss_filter_input.ErrorEllipseMajor_m = std_major(gnss_valid_selector)';
sim_gnss_filter_input.ErrorEllipseMinor_m = std_minor(gnss_valid_selector)';
sim_gnss_filter_input.ErrorEllipseOrientation_deg = alpha(gnss_valid_selector)';
sim_gnss_filter_input.VelocityGroundSigma_ms = squeeze(simout_sensor_data_P.Data(3,3,gnss_valid_selector));
sim_gnss_filter_input.HeadingSigma_deg = rad2deg(squeeze(simout_sensor_data_P.Data(4,4,gnss_valid_selector)));

sim_imu_filter_input = table();
sim_imu_filter_input.Time = simulink_time;
sim_imu_filter_input.UtcTime = datetime(t0+sim_imu_filter_input.Time,'ConvertFrom','posixtime');
sim_imu_filter_input.AccX_mss = simout_sensor_data.Data(:,5);
sim_imu_filter_input.TurnRateZ_degs = simout_sensor_data.Data(:,6);
sim_imu_filter_input.AccXSigma_mss = squeeze(simout_sensor_data_P.Data(5,5,:));
sim_imu_filter_input.TurnRateZSigma_degs = rad2deg(squeeze(simout_sensor_data_P.Data(6,6,:)));

fprintf('done!\n');

%% EKF Filter Data

fprintf('\t 4/%i ''EKF data''...',num_processing_steps);

% Calculations ____________________________________________________________
x_utm = p_0_utm(1)+simout_ekf_x.Data(:,1);
y_utm = p_0_utm(2)+simout_ekf_x.Data(:,2);
[lat,lon] = utm2ll(x_utm,y_utm,p_0_utm(3),'wgs84');

yaw_deg = simplifyHeadingD(rad2deg(simout_ekf_x.Data(:,6)),'180');
v_vehicle = simout_ekf_x.Data(:,4);
heading_deg = yaw_deg;
heading_deg(v_vehicle<0) = simplifyHeadingD(heading_deg(v_vehicle<0)-180,'180');

major_minor_confidence = chi2cdf(1^2,2); % calculate probability for 1x sigma 2D-error-ellipse
[std_major,std_minor,alpha] = calcErrorEllipseParameters(simout_ekf_P.Data(1:2,1:2,:),major_minor_confidence);

P_pos = reshape(simout_ekf_P.Data(1:2,1:2,:),2,2*length(simout_ekf_P.Time));
P_pos = mat2cell(P_pos,2,repmat(2,1,length(simout_ekf_P.Time)))';

d_vehicle = cumsum(abs([0;diff(simout_ekf_x.Data(:,3))]));

% Write to table __________________________________________________________
sim_ekf = table();
sim_ekf.Time = simulink_time;
sim_ekf.UtcTime = datetime(t0+sim_ekf.Time,'ConvertFrom','posixtime');
sim_ekf.UtmEast_m = x_utm;
sim_ekf.UtmNorth_m = y_utm;
sim_ekf.Latitude_deg = lat;
sim_ekf.Longitude_deg = lon;
sim_ekf.DistanceVehicle_m = d_vehicle;
sim_ekf.DistanceTrack_m = simout_ekf_x.Data(:,3);
sim_ekf.VelocityVehicle_ms = v_vehicle;
sim_ekf.AccelerationVehicle_ms = simout_ekf_x.Data(:,5);
sim_ekf.Heading_deg = heading_deg;
sim_ekf.Yaw_deg = yaw_deg;
sim_ekf.TurnRateZ_degs = simout_ekf_x.Data(:,7);
sim_ekf.LatitudeSigma_m = squeeze(simout_ekf_P.Data(2,2,:));
sim_ekf.LongitudeSigma_m = squeeze(simout_ekf_P.Data(1,1,:));
sim_ekf.PositionCov = P_pos;
sim_ekf.ErrorEllipseMajor_m = std_major(:);
sim_ekf.ErrorEllipseMinor_m = std_minor(:);
sim_ekf.ErrorEllipseOrientation_deg = alpha(:);
sim_ekf.DistanceTrackSigma_m = squeeze(simout_ekf_P.Data(3,3,:));
sim_ekf.VelocityVehicleSigma_ms = squeeze(simout_ekf_P.Data(4,4,:));
sim_ekf.AccelerationVehicleSigma_ms = squeeze(simout_ekf_P.Data(5,5,:));
sim_ekf.HeadingSigma_deg = rad2deg(squeeze(simout_ekf_P.Data(6,6,:)));
sim_ekf.TurnRateZSigma_degs = rad2deg(squeeze(simout_ekf_P.Data(7,7,:)));

fprintf('done!\n');

%% IMM Filter Data

fprintf('\t 4/%i ''IMM data''...',num_processing_steps);

% Calculations ____________________________________________________________
x_utm = p_0_utm(1)+simout_imm_x.Data(:,1);
y_utm = p_0_utm(2)+simout_imm_x.Data(:,2);
[lat,lon] = utm2ll(x_utm,y_utm,p_0_utm(3),'wgs84');

yaw_deg = simplifyHeadingD(rad2deg(simout_imm_x.Data(:,6)),'180');
v_vehicle = simout_imm_x.Data(:,4);
heading_deg = yaw_deg;
heading_deg(v_vehicle<0) = simplifyHeadingD(heading_deg(v_vehicle<0)-180,'180');

major_minor_confidence = chi2cdf(1^2,2); % calculate probability for 1x sigma 2D-error-ellipse 
[std_major,std_minor,alpha] = calcErrorEllipseParameters(simout_imm_P.Data(1:2,1:2,:),major_minor_confidence);

P_pos = reshape(simout_imm_P.Data(1:2,1:2,:),2,2*length(simout_imm_P.Time));
P_pos = mat2cell(P_pos,2,repmat(2,1,length(simout_imm_P.Time)))';

d_vehicle = cumsum(abs([0;diff(simout_imm_x.Data(:,3))]));

% Write to table __________________________________________________________
sim_imm = table();
sim_imm.Time = simulink_time;
sim_imm.UtcTime = datetime(t0+sim_imm.Time,'ConvertFrom','posixtime');
sim_imm.UtmEast_m = x_utm;
sim_imm.UtmNorth_m = y_utm;
sim_imm.Latitude_deg = lat;
sim_imm.Longitude_deg = lon;
sim_imm.DistanceVehicle_m = d_vehicle;
sim_imm.DistanceTrack_m = simout_imm_x.Data(:,3);
sim_imm.VelocityVehicle_ms = v_vehicle;
sim_imm.AccelerationVehicle_ms = simout_imm_x.Data(:,5);
sim_imm.Heading_deg = heading_deg;
sim_imm.Yaw_deg = yaw_deg;
sim_imm.TurnRateZ_degs = simout_imm_x.Data(:,7);
sim_imm.LatitudeSigma_m = squeeze(simout_imm_P.Data(2,2,:));
sim_imm.LongitudeSigma_m = squeeze(simout_imm_P.Data(1,1,:));
sim_imm.PositionCov = P_pos;
sim_imm.ErrorEllipseMajor_m = std_major(:);
sim_imm.ErrorEllipseMinor_m = std_minor(:);
sim_imm.ErrorEllipseOrientation_deg = alpha(:);
sim_imm.DistanceTrackSigma_m = squeeze(simout_imm_P.Data(3,3,:));
sim_imm.VelocityVehicleSigma_ms = squeeze(simout_imm_P.Data(4,4,:));
sim_imm.AccelerationVehicleSigma_ms = squeeze(simout_imm_P.Data(5,5,:));
sim_imm.HeadingSigma_deg = rad2deg(squeeze(simout_imm_P.Data(6,6,:)));
sim_imm.TurnRateZSigma_degs = rad2deg(squeeze(simout_imm_P.Data(7,7,:)));
sim_imm.ModelProbabilities = simout_imm_mu.Data;

fprintf('done!\n');

%% TGC Data

fprintf('\t 5/%i ''TGC data''...',num_processing_steps);

% Calculations ____________________________________________________________
x_utm = p_0_utm(1)+simout_loc_position.Data(:,1);
y_utm = p_0_utm(2)+simout_loc_position.Data(:,2);
[lat,lon] = utm2ll(x_utm,y_utm,p_0_utm(3),'wgs84');

yaw_deg = simplifyHeadingD(rad2deg(simout_loc_position.Data(:,4)),'180');
% yaw_deg = rad2deg(simout_loc_position.Data(:,4));
v_vehicle = simout_loc_position.Data(:,3);
heading_deg = yaw_deg;
heading_deg(v_vehicle<0) = simplifyHeadingD(heading_deg(v_vehicle<0)-180,'180');
% heading_deg(v_vehicle<0) = heading_deg(v_vehicle<0)-180;

major_minor_confidence = chi2cdf(1^2,2); % calculate probability for 1x sigma 2D-error-ellipse 
[std_major,std_minor,alpha] = calcErrorEllipseParameters(simout_loc_P.Data(1:2,1:2,:),major_minor_confidence);

P_pos = reshape(simout_loc_P.Data(1:2,1:2,:),2,2*length(simout_loc_P.Time));
P_pos = mat2cell(P_pos,2,repmat(2,1,length(simout_loc_P.Time)))';

% Write to table __________________________________________________________
sim_tgc = table();
sim_tgc.Time = simulink_time;
sim_tgc.UtcTime = datetime(t0+sim_tgc.Time,'ConvertFrom','posixtime');
sim_tgc.UtmEast_m = x_utm;
sim_tgc.UtmNorth_m = y_utm;
sim_tgc.Latitude_deg = lat;
sim_tgc.Longitude_deg = lon;
sim_tgc.DistanceVehicle_m = sim_imm.DistanceVehicle_m;
sim_tgc.DistanceTrack_m = sim_imm.DistanceTrack_m;
sim_tgc.VelocityVehicle_ms = v_vehicle;
sim_tgc.Heading_deg = heading_deg;
sim_tgc.Yaw_deg = yaw_deg;
sim_tgc.LatitudeSigma_m = squeeze(simout_loc_P.Data(2,2,:));
sim_tgc.LongitudeSigma_m = squeeze(simout_loc_P.Data(1,1,:));
sim_tgc.PositionCov = P_pos;
sim_tgc.ErrorEllipseMajor_m = std_major(:);
sim_tgc.ErrorEllipseMinor_m = std_minor(:);
sim_tgc.ErrorEllipseOrientation_deg = alpha(:);
sim_tgc.VelocityVehicleSigma_ms = squeeze(simout_loc_P.Data(3,3,:));
sim_tgc.HeadingSigma_deg = rad2deg(squeeze(simout_loc_P.Data(4,4,:)));

fprintf('done!\n');

%% INS Data

fprintf('\t 6/%i ''INS data''...',num_processing_steps);

% Reference speeds
gnss_v_ground = sqrt(sim_gnss.VelocityNorth_ms.^2+sim_gnss.VelocityEast_ms.^2);
ins_v_ground = sqrt(sim_v_ned.VelocityNorth_ms.^2+sim_v_ned.VelocityEast_ms.^2);

% Reset speed to zero where standstill has been detected
cs_acc_uncorrected = sim_a_ib_b.AccX_mss;% -sim_a_bias.BiasAccX_mss;
cs_acc_uncorrected(logical(sim_standstill_flag.standstill_flag)) = 0; 

% cs_acc_x_corrected = sim_a_eb_n.AccX_mss;
% cs_acc_x_corrected(logical(sim_standstill_flag.standstill_flag)) = 0;      
% cs_acc_y_corrected = sim_a_eb_n.AccY_mss;
% cs_acc_y_corrected(logical(sim_standstill_flag.standstill_flag)) = 0; 

cs_acc_x_corrected = sim_a_eb_n.AccX_mss;
cs_acc_y_corrected = sim_a_eb_n.AccY_mss;
cs_acc_corrected = sqrt(cs_acc_x_corrected.^2+cs_acc_y_corrected.^2);

psi_ins = atan2d(cs_acc_y_corrected,cs_acc_x_corrected);
psi_heading = sim_attitude.Heading_deg;
delta_psi = abs(angdiff(deg2rad(psi_ins),deg2rad(psi_heading)));
dir_neg_selector = (delta_psi > pi/2);

cs_acc_x_corrected(logical(sim_standstill_flag.standstill_flag)) = 0;
cs_acc_y_corrected(logical(sim_standstill_flag.standstill_flag)) = 0;
cs_acc_corrected(dir_neg_selector) = -1*cs_acc_corrected(dir_neg_selector);
cs_acc_corrected(logical(sim_standstill_flag.standstill_flag)) = 0;
    
sim_cs_v_ground_uncorrected = zeros(length(sim_a_ib_b.Time),1);
sim_cs_v_ground_x_corrected = zeros(length(sim_a_eb_n.Time),1);
sim_cs_v_ground_y_corrected = zeros(length(sim_a_eb_n.Time),1);
sim_cs_v_ground_corrected = zeros(length(sim_a_eb_n.Time),1);
standstill_idx = strfind(sprintf('%d',sim_standstill_flag.standstill_flag),'01');
standstill_idx = [standstill_idx,length(sim_standstill_flag.Time)];
for i = 1:length(standstill_idx)  
    
    if i > 1
        drive_idx = standstill_idx(i-1)+1:standstill_idx(i);
    else
        drive_idx = 1:standstill_idx(i);
    end % if
    
    if length(drive_idx)==1
        continue
    end % if
    
    % Manual integration of raw acc data      
    sim_cs_v_ground_uncorrected(drive_idx) = cumtrapz(sim_a_ib_b.Time(drive_idx)-sim_a_ib_b.Time(drive_idx(1)),cs_acc_uncorrected(drive_idx));

    % Manual integration strapdown corrected acc data
    sim_cs_v_ground_x_corrected(drive_idx) = cumtrapz(sim_a_eb_n.Time(drive_idx)-sim_a_eb_n.Time(drive_idx(1)),cs_acc_x_corrected(drive_idx));
    sim_cs_v_ground_y_corrected(drive_idx) = cumtrapz(sim_a_eb_n.Time(drive_idx)-sim_a_eb_n.Time(drive_idx(1)),cs_acc_y_corrected(drive_idx));
    sim_cs_v_ground_corrected(drive_idx) = cumtrapz(sim_a_eb_n.Time(drive_idx)-sim_a_eb_n.Time(drive_idx(1)),cs_acc_corrected(drive_idx));
     
%     cs_v_ground_x_corrected = cumtrapz(sim_a_eb_n.Time(drive_idx),cs_acc_x_corrected(drive_idx));           
%     cs_v_ground_y_corrected = cumtrapz(sim_a_eb_n.Time(drive_idx),cs_acc_y_corrected(drive_idx));
%     sim_cs_v_ground_corrected(drive_idx) = sqrt(cs_v_ground_x_corrected.^2+cs_v_ground_y_corrected.^2);    
    
%     sim_cs_v_ground_uncorrected(i+1:end) = sim_cs_v_ground_uncorrected(i+1:end)-sim_cs_v_ground_uncorrected(i+1);
%     sim_cs_v_ground_corrected(i+1:end) = sim_cs_v_ground_corrected(i+1:end)-sim_cs_v_ground_corrected(i+1);
end % for i
sim_cs_v_ground_xy_corrected = sqrt(sim_cs_v_ground_x_corrected.^2+sim_cs_v_ground_y_corrected.^2);

fprintf('done!\n');

%% Error Data

fprintf('\t 7/%i ''Error data''...',num_processing_steps);

calcErrorData

fprintf('done!\n');

%% Railway-map Data

fprintf('\t 8/%i ''Railway-map data''...\n',num_processing_steps);

% Calculations ____________________________________________________________
%num_tracks = sum(~isnan(simout_topology.signals.values(:,1,end)));
num_tracks = sum(~isnan(simout_track_start_points.signals.values(:,1)));
num_track_elements = max(sum(~isnan(simout_track_maps.signals.values(:,1:end,end))));

% sim_topology = simout_topology.signals.values(1:num_tracks,1:num_tracks,end);
sim_topology = [ zeros(num_tracks-1,1),     eye(num_tracks-1); ... 
                            zeros(1,1), zeros(1,num_tracks-1)];
sim_track_start_points = simout_track_start_points.signals.values(1:num_tracks,1:end,end);
sim_track_maps = simout_track_maps.signals.values(1:num_track_elements,1:end,end);

% sim_track_start_points_cov = cell(num_tracks,1);
% for i = 1:num_tracks
%     sim_track_start_points_cov{i} = simout_P_track_start_points.signals.values(:,:,i,end);
% end % for i
sim_track_start_points_cov = reshape(simout_P_track_start_points.signals.values(:,:,1:num_tracks,end),3,3*num_tracks);
sim_track_start_points_cov = mat2cell(sim_track_start_points_cov,3,repmat(3,1,num_tracks))';

% sim_track_maps_cov = cell(num_track_elements,1);
% for i = 1:num_track_elements
%     sim_track_maps_cov{i} = simout_P_track_maps.signals.values(:,:,i,end);
% end % for i
sim_track_maps_cov = reshape(simout_P_track_maps.signals.values(:,:,1:num_track_elements,end),3,3*num_tracks);
sim_track_maps_cov = mat2cell(sim_track_maps_cov,3,repmat(3,1,num_track_elements))';

sim_map.topology = sim_topology;
sim_map.track_start_points = matTrackStartPoints2tableTrackStartPoints(sim_track_start_points,sim_track_start_points_cov);
sim_map.track_start_points.phi_0 = simplifyHeadingD(sim_map.track_start_points.phi_0,'180');
% sim_map.track_start_points.phi_0 = mod(sim_map.track_start_points.phi_0,360);
sim_map.track_maps = matTrackMap2tableTrackMap(sim_track_maps,sim_track_maps_cov);

fprintf('\t\t')
[sim_map_track_id,sim_map_track_rel_position,sim_map_abs_position_utm,sim_map_orientation,sim_map_curvature,sim_map_radius,sim_map_speed_limit,sim_map_pdated_railway_map] ... 
    = calcMapProperties(sim_map,sim_map_density);
% [sim_map_track_id,sim_map_track_rel_position,sim_map_abs_position_utm,sim_map_orientation,sim_map_curvature,sim_map_radius,sim_map_speed_limit] ... 
%     = calcMapPropertiesFast(sim_map,sim_map_density);
[sim_map_abs_position_lla_latitude,sim_map_abs_position_lla_longitude] = ... 
    utm2ll(p_0_utm(1)+sim_map_abs_position_utm(1,:),p_0_utm(2)+sim_map_abs_position_utm(2,:),p_0_utm(3),'wgs84');

straight_tm_selector = (sim_map.track_maps.track_element == 1);
straight_sp_selector = ismember(sim_map.track_start_points.ID,unique(sim_map.track_maps(straight_tm_selector,:).ID));
sim_map_straights.topology = sim_map.topology(straight_sp_selector,straight_sp_selector);
sim_map_straights.track_start_points = sim_map.track_start_points(straight_sp_selector,:);
sim_map_straights.track_maps = sim_map.track_maps(straight_tm_selector,:);

if ~isempty(sim_map_straights.topology)
    fprintf('\t\t')
    [~,~,sim_map_straights_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_straights,sim_map_density);
    % [~,~,sim_map_straights_abs_position_utm,~,~,~,~] = calcMapPropertiesFast(sim_map_straights,sim_map_density);
    [sim_map_straights_abs_position_lla_latitude,sim_map_straights_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+sim_map_straights_abs_position_utm(1,:),p_0_utm(2)+sim_map_straights_abs_position_utm(2,:),p_0_utm(3),'wgs84');
end % if

arc_tm_selector = (sim_map.track_maps.track_element == 3);
arc_sp_selector = ismember(sim_map.track_start_points.ID,unique(sim_map.track_maps(arc_tm_selector,:).ID));
sim_map_arcs.topology = sim_map.topology(arc_sp_selector,arc_sp_selector);
sim_map_arcs.track_start_points = sim_map.track_start_points(arc_sp_selector,:);
sim_map_arcs.track_maps = sim_map.track_maps(arc_tm_selector,:);

if ~isempty(sim_map_arcs.topology)
    fprintf('\t\t')
    [~,~,sim_map_arcs_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_arcs,sim_map_density);
    % [~,~,sim_map_arcs_abs_position_utm,~,~,~,~] = calcMapPropertiesFast(sim_map_arcs,sim_map_density);
    [sim_map_arcs_abs_position_lla_latitude,sim_map_arcs_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+sim_map_arcs_abs_position_utm(1,:),p_0_utm(2)+sim_map_arcs_abs_position_utm(2,:),p_0_utm(3),'wgs84');
end % if



% Write to table __________________________________________________________
sim_tgc_map = table();
sim_tgc_map.track_id = sim_map_track_id;
sim_tgc_map.rel_position = sim_map_track_rel_position;
sim_tgc_map.orientation = {sim_map_orientation};
sim_tgc_map.curvature = sim_map_curvature;
sim_tgc_map.radius = sim_map_radius;
sim_tgc_map.speed_limit = sim_map_speed_limit;
sim_tgc_map.UtmEast_m = sim_map_abs_position_utm(1,:)+p_0_utm(1);
sim_tgc_map.UtmNorth_m = sim_map_abs_position_utm(2,:)+p_0_utm(2);
sim_tgc_map.Latitude_deg = sim_map_abs_position_lla_latitude;
sim_tgc_map.Longitude_deg = sim_map_abs_position_lla_longitude;  

if ~isempty(sim_map_straights.topology)
    sim_tgc_map.Straights_UtmEast_m = sim_map_straights_abs_position_utm(1,:)+p_0_utm(1);
    sim_tgc_map.Straights_UtmNorth_m = sim_map_straights_abs_position_utm(2,:)+p_0_utm(2);
    sim_tgc_map.Straights_Latitude_deg = sim_map_straights_abs_position_lla_latitude;
    sim_tgc_map.Straights_Longitude_deg = sim_map_straights_abs_position_lla_longitude; 
end % if

if ~isempty(sim_map_arcs.topology)
    sim_tgc_map.Arcs_UtmEast_m = sim_map_arcs_abs_position_utm(1,:)+p_0_utm(1);
    sim_tgc_map.Arcs_UtmNorth_m = sim_map_arcs_abs_position_utm(2,:)+p_0_utm(2);
    sim_tgc_map.Arcs_Latitude_deg = sim_map_arcs_abs_position_lla_latitude;
    sim_tgc_map.Arcs_Longitude_deg = sim_map_arcs_abs_position_lla_longitude; 
end % if

fprintf('\t done!\n');

%% Optimization data

fprintf('\t 9/%i ''Optimization data''...',num_processing_steps);

gnss_valid_selector = logical(simout_gnss_valid_flag.Data(:,1));
standstill_idx = (sim_imm.VelocityVehicle_ms == 0);

% Map-Optimization data (IMM) _____________________________________________

optimization_nonnan_data_selector_tracks = ~isnan(simout_optimization_data_selector.signals.values(:,1,end));
optimization_data_selector_temp = simout_optimization_data_selector.signals.values(optimization_nonnan_data_selector_tracks,:,end);
optimization_imm_data_selector = [optimization_data_selector_temp(:,1), zeros(size(optimization_data_selector_temp,1),optimization_data_selector_temp(end,2))];
for i = 1:size(optimization_data_selector_temp,1)
    if i == 1
        selector_start_index = 2;
    else
        selector_start_index = optimization_data_selector_temp(i-1,2)+2;
    end % if        
    selector_end_index = optimization_data_selector_temp(i,2)+1;        
    optimization_imm_data_selector(i,selector_start_index:selector_end_index) = 1;
end % for i

% Map-Optimization data (GPS) _____________________________________________
optimization_gps_data_selector = optimization_imm_data_selector(:,[true;gnss_valid_selector]');

fprintf('done!\n');
  
%% Finish

simout_data_prepared = true;

%% Save

prompt_str = 'Do you want to save the processed simout data (j/n)?';
user_str = input(prompt_str,'s');
% user_str = 'j';

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
    