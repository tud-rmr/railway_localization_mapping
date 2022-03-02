% 
% Configuration file for localization with the data set 'Chemnitz' (C)
% (see: https://doi.org/10.25534/tudatalib-360)
%   
%   Author: Hanno Winter
%   Date: 08-Apr-2021; Last revision: 14-Dec-2021

%% Data

% Normal Driving (short) __________________________________________________
% start_time = '2018-10-24 11:50:20'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2018-10-24 11:52:30'; % format: uuuu-MM-dd HH:mm:ss

% Normal Driving __________________________________________________________
% start_time = '2018-10-24 12:08:15'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2018-10-24 12:23:00'; % format: uuuu-MM-dd HH:mm:ss

% Reverse Driving _________________________________________________________
% start_time = '2018-10-24 13:19:45'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2018-10-24 13:23:40'; % format: uuuu-MM-dd HH:mm:ss

% Makersbach ______________________________________________________________
% start_time = '2018-10-24 13:43:00'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2018-10-24 13:57:00'; % format: uuuu-MM-dd HH:mm:ss

% Multiple starts and stops with reversing ________________________________
% start_time = '2018-10-24 14:01:30'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2018-10-24 14:27:25'; % format: uuuu-MM-dd HH:mm:ss

% Schwarzenberg to Annaberg-Buchholz Süd (whole section with reference map available)
start_time = '2018-10-24 14:20:00'; % format: uuuu-MM-dd HH:mm:ss
end_time = '2018-10-24 15:24:21'; % format: uuuu-MM-dd HH:mm:ss

%% Simulation 

% Sample times ____________________________________________________________
imu_sample_time = 0.002; % sampe time of imu data in sec
gnss_sample_time = 1; % sample time of gnss data in sec
imm_sample_time = 0.1; % sample time of IMM filter in sec
imm_mu_sample_time = 0.25; % sample time for model update in IMM filter in sec

%% Data preparation

max_pdop = []; % Max GPS measurement PDOP (set empty if all available measurement data should be used)

%% Strapdown

f_euler_ma_time = 0.2;
comp_filter_static_gain = 0.02;
comp_filter_adaptive_gain_thresholds = [0 5e-6];
E_w_limit = 2.5e-6;
E_w_ma_time = 3;
a_e_limit = 0.07;
a_e_ma_time = 2;
w_bias_gain = 1e-3;
a_bias_gain = 1e-3;
gnss_speed_limit = 0.1;

d_bogies = 15.1; % distance between bogies in m (set to '0' if not relevant)

%% IMM

% Model parameters ________________________________________________________
arc_w_min = 0.0035; % Minimum turn-rate for circular-arc model

% Model transition probabilites ___________________________________________
% 
%   Model order: [straight, arc, trash];
% 
Pi = [0.900 0.000 0.100; ... 
      0.000 0.900 0.100; ... 
      0.300 0.100 0.600];
  
% Initial model probabilities _____________________________________________
% 
%   Model order: [straight, arc, trash];
% 
mu_0 = [0.34; 0.33; 0.33];

% Initial uncertainties ___________________________________________________
sigma_d_0 = 0.1; % standard deviation of the initial distance in m
sigma_psi_0 = deg2rad(180); % standard deviation of the initial heading in rad

% IMU measurement noise ___________________________________________________
sigma_v_acc  = 0.1*10; % 9.81;  % standard deviation of acceleration sensor (in m/s^2)
sigma_v_gyro = 0.0020; % deg2rad(0.1);  % standard deviation of turn rate sensor (in rad/s)

% System noise_____________________________________________________________
sigma_w_imm_trash_pos = 7*imm_sample_time; % standard deviation of postion estimate 
sigma_w_imm_trash_acc  = 1.5*imm_sample_time;  % standard deviation of system motion
sigma_w_imm_trash_gyro = 0.0400*imm_sample_time; % standard deviation of system turn rate

sigma_w_imm_linear_pos = 7*imm_sample_time; % standard deviation of postion estimate 
sigma_w_imm_linear_acc  = 1.5*imm_sample_time; % standard deviation of system motion
sigma_w_imm_linear_gyro = 0.0020*imm_sample_time;   % standard deviation of system turn rate

sigma_w_imm_circle_pos = 7*imm_sample_time; % standard deviation of postion estimate 
sigma_w_imm_circle_acc  = 1.5*imm_sample_time;  % standard deviation of system motion
sigma_w_imm_circle_gyro = 0.0020*imm_sample_time; % standard deviation of system turn rate

%% TGC

% Min speed _______________________________________________________________
v_min_tgc_speed = round(2*1/3.6 * 10e2)/10e2; % Minimum speed to perform track geometry position correction
v_min_tgc_hysteresis = round(0.4/3.6 * 10e2)/10e2; % Hysteresis for minimum TGC speed enabling/disabling

% Variables sizes _________________________________________________________
max_track_map_size = 500; % max size of track-maps (needed for code generation within Simulink)

% System uncertainties during parameter estimation ________________________
sigma_psi_straight = deg2rad(4)*imm_sample_time;
sigma_psi_arc = deg2rad(4)*imm_sample_time;
sigma_p_0 = 0.1*imm_sample_time;
sigma_r = 0.1*imm_sample_time;

%% Set parameters

% Make sure that parameters are updated, by resetting init bit (init 
% script should be started manually or by the simulink model)
init_localization_executed = false;
