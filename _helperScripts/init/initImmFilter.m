% 
% Initialize IMM module
%   
%   Author: Hanno Winter
%   Date: 20-Mar-2021; Last revision: 17-Apr-2021

%% Init

if ~exist('init_localization_executed','var')
    warning('initImmFilter.m: This script is supposed to be called from ''initLocalization.m''!')
    % return
elseif ~exist('prepare_sim_input_data_executed','var')
    warning('initImmFilter.m: Run ''prepareSimInputData.m'' before executing this script!')
    % return
else
    init_c_strapdown_executed = false;
end % if

%% Settings

% System and measurement functions ________________________________________
imm_state_transition_fcns = { 'simFcn_Imm_Linear_StateTransitionFcn', ... 
                              'simFcn_Imm_Circle_StateTransitionFcn', ... 
                              'simFcn_Imm_Trash_StateTransitionFcn'};
imm_measurement_fcns = {'simFcn_Imm_MeasurementFcn', ... 
                        'simFcn_Imm_MeasurementFcn', ... 
                        'simFcn_Imm_MeasurementFcn'};
imm_residuum_fcns = {'simFcn_Imm_ResiduumFcn', ... 
                     'simFcn_Imm_ResiduumFcn', ... 
                     'simFcn_Imm_ResiduumFcn'};

%% Initial states

p_0   = zeros(2,1); % in m
d_0   = 0; % in m
v_0   = sqrt(sum(v_ned_eb_n_0([1,2]).^2)); % in m/s
a_0   = 0; % sqrt(sum(acc_data.Data(1,[1,2]).^2));  % in m/s^2
psi_0 = deg2rad(heading_0); % in rad
w_0   = 0; % deg2rad(gyro_data.Data(1,3)); % in rad/s

x_0   = [p_0;d_0;v_0;a_0;psi_0;w_0];
x_0_imm = repmat(x_0,1,3);

%% Initial covariance matrices

% sigma_w_imm_linear = blkdiag(sigma_w_imm_linear_acc,sigma_w_imm_linear_gyro);
% sigma_w_imm_circle = blkdiag(sigma_w_imm_circle_acc,sigma_w_imm_circle_gyro);
% sigma_w_imm_trash = blkdiag(sigma_w_imm_trash_acc,sigma_w_imm_trash_gyro);
% 
sigma_w_imm_linear = blkdiag(sigma_w_imm_linear_pos,sigma_w_imm_linear_pos,sigma_w_imm_linear_acc,sigma_w_imm_linear_gyro);
sigma_w_imm_circle = blkdiag(sigma_w_imm_circle_pos,sigma_w_imm_circle_pos,sigma_w_imm_circle_acc,sigma_w_imm_circle_gyro);
sigma_w_imm_trash = blkdiag(sigma_w_imm_trash_pos,sigma_w_imm_trash_pos,sigma_w_imm_trash_acc,sigma_w_imm_trash_gyro);
% 
% Gamma_linear = [[1/2*imm_sample_time^2 1/2*imm_sample_time^2 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', [0 0 0 0 0 imm_sample_time 1]'];
% Gamma_circle = [[1/2*imm_sample_time^2 1/2*imm_sample_time^2 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', [0 0 0 0 0 imm_sample_time 1]'];
% Gamma_trash  = [[1/2*imm_sample_time^2 1/2*imm_sample_time^2 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', [0 0 0 0 0 imm_sample_time 1]'];
% 
% Gamma_linear = [[0 0 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', [0 0 0 0 0 imm_sample_time 1]'];
% Gamma_circle = [[0 0 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', [0 0 0 0 0 imm_sample_time 1]'];
% Gamma_trash  = [[0 0 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', [0 0 0 0 0 imm_sample_time 1]'];
% 
Gamma_linear = [ ... 
                 [1 0 0 0 0 0 0]', ... 
                 [0 1 0 0 0 0 0]', ... 
                 [0 0 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', ... 
                 [0 0 0 0 0 imm_sample_time 1]' ... 
               ];
Gamma_circle = [ ... 
                 [1 0 0 0 0 0 0]', ... 
                 [0 1 0 0 0 0 0]', ... 
                 [0 0 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', ... 
                 [0 0 0 0 0 imm_sample_time 1]' ... 
               ];
Gamma_trash = [ ... 
                 [1 0 0 0 0 0 0]', ... 
                 [0 1 0 0 0 0 0]', ... 
                 [0 0 1/2*imm_sample_time^2 imm_sample_time 1 0 0]', ... 
                 [0 0 0 0 0 imm_sample_time 1]' ... 
               ];


Q_imm_linear = Gamma_linear*sigma_w_imm_linear.^2*Gamma_linear';
% Q_imm_linear(1:2,1:2) = eye(2)*sigma_w_imm_linear_pos^2;

Q_imm_circle = Gamma_circle*sigma_w_imm_circle.^2*Gamma_circle';
% Q_imm_circle(1:2,1:2) = eye(2)*sigma_w_imm_circle_pos^2;

Q_imm_trash = Gamma_trash*sigma_w_imm_trash.^2*Gamma_trash';
% Q_imm_trash(1:2,1:2) = eye(2)*sigma_w_imm_trash_pos^2;

% Q_imm_linear = floor(Q_imm_linear.*1e8)/1e8;
% Q_imm_circle = floor(Q_imm_circle.*1e8)/1e8;
% Q_imm_trash = floor(Q_imm_trash.*1e8)/1e8;

% Q_imm_linear = diag(([0.1 0.1 0.1 1 1.5 0.0002 0.002]*imm_sample_time).^2);
% Q_imm_circle = diag(([0.1 0.1 0.1 1 1.5 0.0002 0.002]*imm_sample_time).^2);
% Q_imm_trash  = diag(([0.1 0.1 0.1 1 1.5 0.0040 0.040]*imm_sample_time).^2);

% Q_imm_linear = diag(([10 10 1 1 100 deg2rad(0.1) 0.002]*imm_sample_time).^2);
% Q_imm_circle = diag(([10 10 1 1 100 deg2rad(0.1) 0.002]*imm_sample_time).^2);
% Q_imm_trash  = diag(([10 10 1 1 100 deg2rad(0.1) 0.010]*imm_sample_time).^2);

% Q_imm_linear = diag(([0.01 0.01 0.01 1 1 deg2rad(0.1) 0.002]*imm_sample_time).^2);
% Q_imm_circle = diag(([0.01 0.01 0.01 1 1 deg2rad(0.1) 0.002]*imm_sample_time).^2);
% Q_imm_trash  = diag(([1 1 1 1 1 deg2rad(0.1) 0.010]*imm_sample_time).^2);

Q_imm = cat(3,Q_imm_linear,Q_imm_circle,Q_imm_trash);

P_0 = blkdiag( ... 
               gnss_data.LatitudeSigma_m(1)^2,       ... % positon
               gnss_data.LongitudeSigma_m(1)^2,      ... % position
               sigma_d_0^2,                           ... % distance
               gnss_data.VelocityGroundSigma_ms(1)^2, ... % speed
               sigma_v_acc^2,                        ... % acceleration
               sigma_psi_0^2,                        ... % psi           % deg2rad(gnss_data.HeadingSigma_deg(heading_ok_index))^2
               sigma_v_gyro^2                        ... % omega
             );
P_0_imm = repmat(P_0,1,1,3);
