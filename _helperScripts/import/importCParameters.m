% 
%   Data available at: http://dx.doi.org/10.25534/tudatalib-166.2
%
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 18-Apr-2021
% 

%% Init 

loadCFilenamePatterns

%% Checks

% if ~exist('import_C_sensor_parameters_executed','var')
%     import_C_sensor_parameters_executed = false;
% elseif import_C_sensor_parameters_executed
%     fprintf('Import parameters script already executed!\n');
%     return  
% end % if

%% Parameter files
% 
%   Data available at: http://dx.doi.org/10.25534/tudatalib-166.2
%                       

gnss_thales_parameters_paths = { ... 
                                 fullfile('Dataset_C_20181024','_parameters','gnss_thales_parameters.csv') ... 
                               };

gnss_inatm200stn_parameters_paths = { ... 
                                      fullfile('Dataset_C_20181024','_parameters','gnss_inatm200stn_parameters.csv') ... 
                                    };
                         
imu_inatm200stn_parameters_paths = { ... 
                                     fullfile('Dataset_C_20181024','_parameters','imu_inatm200stn_parameters.csv') ... 
                                   };

special_var_types = { ... 
                      'LeverArmX_m','double'; ... 
                      'LeverArmY_m','double'; ... 
                      'LeverArmZ_m','double'; ... 
                      'MountingRoll_deg','double'; ... 
                      'MountingPitch_deg','double'; ... 
                      'MountingYaw_deg','double' ...
                    };                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Thales (GNSS) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(gnss_thales_parameters_paths,gnss_thales_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,gnss_thales_parameters_root_string,gnss_thales_parameters_root_string);
% end % if
assignin('base',gnss_thales_parameters_root_string,out_data); 
clear out_data
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iNat M200 STN (GNSS) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(gnss_inatm200stn_parameters_paths,gnss_inatm200stn_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,gnss_inatm200stn_parameters_root_string,gnss_inatm200stn_parameters_root_string);
% end % if
assignin('base',gnss_inatm200stn_parameters_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iNat M200 STN (IMU) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(imu_inatm200stn_parameters_paths,imu_inatm200stn_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,imu_inatm200stn_parameters_root_string,imu_inatm200stn_parameters_root_string);
% end % if
assignin('base',imu_inatm200stn_parameters_root_string,out_data); 
clear out_data

%% Finish script

import_C_sensor_parameters_executed = true;
