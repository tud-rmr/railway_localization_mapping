% 
% Initialize strapdown module
%   
%   Author: Hanno Winter
%   Date: 20-Mar-2021; Last revision: 17-Apr-2021

%% Init

if ~exist('init_localization_executed','var')
    warning('initStrapdown.m: This script is supposed to be called from ''initLocalization.m''!')
    %return
elseif ~exist('prepare_sim_input_data_executed','var')
    warning('initStrapdown.m: Run ''prepareSimInputData.m'' before executing this script!')
    %return
else
    init_c_strapdown_executed = false;
end % if

%% Initial values

% Initial attitude
heading_ok_index = find((abs(gnss_data.VelocityEast_ms) > 1.0) & (abs(gnss_data.VelocityNorth_ms) > 1.0),1);
heading_0 = simplifyHeadingD(gnss_data.Heading_deg(heading_ok_index),'180');
% heading_0 = atan2d(gnss_data.VelocityEast_ms(heading_ok_index),gnss_data.VelocityNorth_ms(heading_ok_index));
% heading_0 = 0;
q_b_n_0 = quatFromVect(acc_data.Data(1,:)',[0;0;-1]); % find initial pitch and roll
q_b_n_0 = multQuat(eulerToQuat([0;0;heading_0*pi/180]),q_b_n_0); % set initial heading
% q_b_n_0 = eulerToQuat([-10;-10;heading_0]*pi/180);

% Initial speed
v_ned_eb_n_0 = gnss_data{1,{'VelocityNorth_ms','VelocityEast_ms','VelocityDown_ms'}}';

% Initial ellipsoid parameteres
p_llh_n_0_rad = gnss_data{1,{'Latitude_deg','Longitude_deg','AltitudeEllipsoid_m'}}' .* [pi/180;pi/180;1];
[~,~,Rn_0,Re_0,~] = getWgs84Parameter(p_llh_n_0_rad);

% Initial GNSS lever arm correction
gnss_lever_arm_r = transformVect(gnss_lever_arm,q_b_n_0);
delta_gnss_lever_arm_r = [ ...
                           atan( gnss_lever_arm_r(1)/(Rn_0-gnss_lever_arm_r(3)+p_llh_n_0_rad(3)) ); ... 
                           atan( gnss_lever_arm_r(2)/(Re_0-gnss_lever_arm_r(3)+p_llh_n_0_rad(3)) ); ... 
                           -gnss_lever_arm_r(3) ...
                         ];
                     
% Initial position (corrected by lever arm offset)
p_llh_n_0_rad = p_llh_n_0_rad - delta_gnss_lever_arm_r;
p_llh_n_0_deg = [p_llh_n_0_rad(1:2)*180/pi;p_llh_n_0_rad(3)];

%% Finish

init_c_strapdown_executed = true;
