%
%   Data available at: http://dx.doi.org/10.25534/tudatalib-360
%
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 02-Dec-2020
% 

%% Init 

loadBsFilenamePatterns

%% Checks

% if ~exist('import_BS_sensor_parameters_executed','var')
%     import_BS_sensor_parameters_executed = false;
% elseif import_BS_sensor_parameters_executed
%     fprintf('Import parameters script already executed!\n');
%     return  
% end % if

%% Parameter files
%
%   Data available at: http://dx.doi.org/10.25534/tudatalib-360
%

gnss_javad_parameters_paths = { ... 
                                fullfile('Dataset_BS_20190222','_parameters','parameters_javad.csv') ... 
                              };
                       
imu_xsens_parameters_paths = { ... 
                                fullfile('Dataset_BS_20190222','_parameters','parameters_xsens.csv') ... 
                              };

speedDist_odometer_parameters_paths = { ... 
                                        fullfile('Dataset_BS_20190222','_parameters','parameters_odometer.csv') ... 
                                      };

speedDist_correvit_parameters_paths = { ... 
                                        fullfile('Dataset_BS_20190222','_parameters','parameters_correvit.csv') ... 
                                      };
                                 
speedDist_siemens_parameters_paths = { ... 
                                       fullfile('Dataset_BS_20190222','_parameters','parameters_siemens.csv') ... 
                                     };                              
                                
gnss_inatm200stn_parameters_paths = { ... 
                                      fullfile('Dataset_BS_20190222','_parameters','parameters_inatm200stn_gnss.csv') ... 
                                    };
                         
imu_inatm200stn_parameters_paths = { ... 
                                     fullfile('Dataset_BS_20190222','_parameters','parameters_inatm200stn_imu.csv') ... 
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
%% GNSS: Javad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(gnss_javad_parameters_paths,gnss_javad_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
% 	writeToMatFile(out_data,gnss_javad_parameters_root_string,gnss_javad_parameters_root_string);
% end % if 
assignin('base',gnss_javad_parameters_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GNSS: Xsens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(imu_xsens_parameters_paths,imu_xsens_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,imu_xsens_parameters_root_string,imu_xsens_parameters_root_string);
% end % if
assignin('base',imu_xsens_parameters_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed/Distance: Odometer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(speedDist_odometer_parameters_paths,speedDist_odometer_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,speedDist_odometer_parameters_root_string,speedDist_odometer_parameters_root_string);
% end % if
assignin('base',speedDist_odometer_parameters_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed/Distance: Correvit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(speedDist_correvit_parameters_paths,speedDist_correvit_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,speedDist_correvit_parameters_root_string,speedDist_correvit_parameters_root_string);
% end % if
assignin('base',speedDist_correvit_parameters_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Speed/Distance: Siemens Radar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, ~, load_flag] = importFromFile(speedDist_siemens_parameters_paths,speedDist_siemens_parameters_root_string,'DataType','parameters','SpecialVarTypes',special_var_types);
% if load_flag == 1 % save to .mat for faster access in the future
%     writeToMatFile(out_data,speedDist_siemens_parameters_root_string,speedDist_siemens_parameters_root_string);
% end % if
assignin('base',speedDist_siemens_parameters_root_string,out_data); 
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

import_BS_sensor_parameters_executed = true;
