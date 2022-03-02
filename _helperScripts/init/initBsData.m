% 
% Initialize data from data set 'Chemnitz' (see: https://tudatalib.ulb.tu-darmstadt.de/handle/tudatalib/2530)
%   
%   Author: Hanno Winter
%   Date: 10-Apr-2021; Last revision: 10-Apr-2021

%% Load data

loadBsFilenamePatterns
importBsParameters      
importBsLimitedProcessedData

imu_data = eval(imu_inatm200stn_processed_data_root_string);
imu_parameters = eval(imu_inatm200stn_parameters_root_string);
gnss_data = eval(gnss_inatm200stn_processed_data_root_string);
gnss_parameters = eval(gnss_inatm200stn_parameters_root_string);
