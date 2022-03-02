% 
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 18-Apr-2021
% 

%% Checks

% if ~exist('loadCFilenamePatterns_loaded','var')
%     loadCFilenamePatterns_loaded = false;
% elseif loadCFilenamePatterns_loaded
%     % fprintf('Already loaded filname patterns!\n');
%     return  
% end % if

%% Parameters

gnss_thales_parameters_root_string = 'gnss_thales_parameters';
gnss_inatm200stn_parameters_root_string = 'gnss_inatm200stn_parameters';
imu_inatm200stn_parameters_root_string = 'imu_inatm200stn_parameters';

%% Raw Data

gnss_thales_raw_data_root_string = 'gnss_thales_raw_data';
gnss_inatm200stn_raw_data_root_string = 'gnss_inatm200stn_raw_data';
imu_inatm200stn_raw_data_root_string = 'imu_inatm200stn_raw_data';
ref_inatm200stn_raw_data_root_string = 'ref_inatm200stn_raw_data';

%% Processed Data

gnss_thales_processed_data_root_string = 'gnss_thales_processed_data';
gnss_inatm200stn_processed_data_root_string = 'gnss_inatm200stn_processed_data';
imu_inatm200stn_processed_data_root_string = 'imu_inatm200stn_processed_data';
ref_inatm200stn_processed_data_root_string = 'ref_inatm200stn_internal_ekf_data';

%% Reference Data

ref_inatm200stn_reference_data_root_string = ref_inatm200stn_processed_data_root_string;
ref_ekfFusion_inatm200stn_reference_data_root_string = 'ref_ekfFusion_inatm200stn_data';

%% Finish

loadCFilenamePatterns_loaded = true;
