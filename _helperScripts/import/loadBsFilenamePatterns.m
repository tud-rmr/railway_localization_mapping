% 
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 18-Apr-2021
% 

%% Init 

% if ~exist('loadBsFilenamePatterns_loaded','var')
%     loadBsFilenamePatterns_loaded = false;
% elseif loadBsFilenamePatterns_loaded
%     % fprintf('Already loaded filname patterns!\n');
%     return  
% end % if

%% Parameters

gnss_javad_parameters_root_string = 'gnss_javad_parameters';
imu_xsens_parameters_root_string = 'imu_xsens_parameters';
speedDist_odometer_parameters_root_string = 'speedDist_odometer_parameters';
speedDist_correvit_parameters_root_string = 'speedDist_correvit_parameters';
speedDist_siemens_parameters_root_string = 'speedDist_siemens_parameters';
gnss_inatm200stn_parameters_root_string = 'gnss_inatm200stn_parameters';
imu_inatm200stn_parameters_root_string = 'imu_inatm200stn_parameters';

%% Raw Data

gnss_javad_raw_data_root_string = 'gnss_javad_raw_data';
imu_xsens_raw_data_root_string = 'imu_xsens_raw_data';
speedDist_odometer_raw_data_root_string = 'speedDist_odometer_raw_data';
speedDist_correvit_raw_data_root_string = 'speedDist_correvit_raw_data';
speedDist_siemens_raw_data_root_string = 'speedDist_siemens_raw_data';
gnss_inatm200stn_raw_data_root_string = 'gnss_inatm200stn_raw_data';
imu_inatm200stn_raw_data_root_string = 'imu_inatm200stn_raw_data';
ref_inatm200stn_raw_data_root_string = 'ref_inatm200stn_raw_data';

%% Processed Data

gnss_javad_processed_data_root_string = 'gnss_javad_processed_data';
imu_xsens_processed_data_root_string = 'imu_xsens_processed_data';
speedDist_odometer_processed_data_root_string = 'speedDist_odometer_processed_data';
speedDist_correvit_processed_data_root_string = 'speedDist_correvit_processed_data';
speedDist_siemens_processed_data_root_string = 'speedDist_siemens_processed_data';
gnss_inatm200stn_processed_data_root_string = 'gnss_inatm200stn_processed_data';
imu_inatm200stn_processed_data_root_string = 'imu_inatm200stn_processed_data';
ref_inatm200stn_processed_data_root_string = 'ref_inatm200stn_internal_ekf_data';

%% Reference Data

ref_inatm200stn_reference_data_root_string = ref_inatm200stn_processed_data_root_string;
ref_ekfFusion_javad_xsens_reference_data_root_string = 'ref_ekfFusion_javad_xsens_data';
ref_ekfFusion_inatm200stn_reference_data_root_string = 'ref_ekfFusion_inatm200stn_data';

%% Finish

% fprintf('Filename patterns loaded!\n');
loadBsFilenamePatterns_loaded = true;
