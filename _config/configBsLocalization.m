% 
% Configuration file for localization with the data set 'Braunschweig' (BS)
% (see: https://doi.org/10.25534/tudatalib-360)
%   
%   Author: Hanno Winter
%   Date: 10-Apr-2021; Last revision: 09-Oct-2021

%% Data

% Drive I _________________________________________________________________
% start_time = '2019-02-22 11:19:26'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2019-02-22 11:26:00'; % format: uuuu-MM-dd HH:mm:ss

% Drive II (backwards) ____________________________________________________
% start_time = '2019-02-22 11:26:00'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2019-02-22 11:31:40'; % format: uuuu-MM-dd HH:mm:ss

% Drive III _______________________________________________________________
% start_time = '2019-02-22 11:31:40'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2019-02-22 11:36:50'; % format: uuuu-MM-dd HH:mm:ss

% Drive IV (backwards) ____________________________________________________
% start_time = '2019-02-22 11:36:50'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2019-02-22 11:41:20'; % format: uuuu-MM-dd HH:mm:ss

% Drive V _________________________________________________________________
% start_time = '2019-02-22 11:41:20'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2019-02-22 11:46:15'; % format: uuuu-MM-dd HH:mm:ss

% Full session ____________________________________________________________
start_time = '2019-02-22 11:19:26'; % format: uuuu-MM-dd HH:mm:ss
end_time = '2019-02-22 11:46:15'; % format: uuuu-MM-dd HH:mm:ss

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
comp_filter_adaptive_gain_thresholds = [0 5e-3];
E_w_limit = 6.00e-5;
E_w_ma_time = 0.5;
a_e_limit = 0.35;
a_e_ma_time = 0.5;
w_bias_gain = 1e-3;
a_bias_gain = 1e-3;
gnss_speed_limit = 0.1;

d_bogies = 0; % distance between bogies in m (set to '0' if not relevant)

%% IMM

% Model parameters ________________________________________________________
arc_w_min = 0.0050; % Minimum turn-rate for circular-arc model

% Model transition probabilites ___________________________________________
% 
%   Model order: [straight, arc, trash];
% 
Pi = [0.800 0.000 0.200; ... 
      0.000 0.800 0.200; ... 
      0.175 0.025 0.800];
Pi = [0.950 0.000 0.050; ... 
      0.000 0.200 0.800; ... 
      0.400 0.400 0.200];
Pi = [0.700 0.000 0.300; ... 
      0.000 0.700 0.300; ... 
      0.290 0.010 0.700];
  
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
sigma_v_gyro = 0.007; % deg2rad(0.1);  % standard deviation of turn rate sensor (in rad/s)

% System noise_____________________________________________________________
sigma_w_imm_trash_pos = 7*imm_sample_time; % standard deviation of postion estimate 
sigma_w_imm_trash_acc  = 1.5*imm_sample_time;  % standard deviation of system motion
sigma_w_imm_trash_gyro = 0.1400*imm_sample_time; % standard deviation of system turn rate

sigma_w_imm_linear_pos = 7*imm_sample_time; % standard deviation of postion estimate 
sigma_w_imm_linear_acc  = 1.5*imm_sample_time; % standard deviation of system motion
sigma_w_imm_linear_gyro = 0.0070*imm_sample_time;   % standard deviation of system turn rate

sigma_w_imm_circle_pos = 7*imm_sample_time; % standard deviation of postion estimate (bigger than for linear model because for this vehicle the turn rates are very noisy)
sigma_w_imm_circle_acc  = 1.5*imm_sample_time;  % standard deviation of system motion
sigma_w_imm_circle_gyro = 0.0070*imm_sample_time; % standard deviation of system turn rate

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
