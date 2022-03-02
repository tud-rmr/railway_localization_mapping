% 
%   Data available at: http://dx.doi.org/10.25534/tudatalib-166.2
%
%   Author: Hanno Winter
%   Date: 02-Dec-2020; Last revision: 18-Apr-2021
% 

%% Init

% Settings ________________________________________________________________

% start_time = '2018-10-24 14:20:30'; % format: uuuu-MM-dd HH:mm:ss
% end_time = '2018-10-24 15:25:30'; % format: uuuu-MM-dd HH:mm:ss

% Load scripts ____________________________________________________________

loadCFilenamePatterns

%% Checks

% if ~exist('import_C_limited_processed_data_executed','var')
%     import_C_limited_processed_data_executed = false;
% elseif import_C_limited_processed_data_executed
%     fprintf('Import limited processed data script already executed!\n');
%     return  
% end % if

%% Processed Data
% 
%   Data available at: http://dx.doi.org/10.25534/tudatalib-166.2
%
%   Filepath or filename has to contain session information, i.e. it has to 
%   contain a string which can be found with the regular 
%   expression '[Ss]ession\d{2}'. Valid strings are, 
%   e.g. 'Session01','session01','session02', etc
%
%   If the data has been chunked into several parts it is import to provide
%   the order of the parts.
%   

gnss_thales_processed_data_paths = { ... 
                                     fullfile('Dataset_C_20181024','2018-10-24_Session01','02_processed','gnss_thales_processed_data_session01.csv'); ... 
                                     fullfile('Dataset_C_20181024','2018-10-24_Session02','02_processed','gnss_thales_processed_data_session02.csv'); ... 
                                     fullfile('Dataset_C_20181024','2018-10-24_Session03','02_processed','gnss_thales_processed_data_session03.csv'); ... 
                                   };

gnss_inatm200stn_processed_data_paths = { ...                                           
                                          fullfile('Dataset_C_20181024','2018-10-24_Session01','02_processed','gnss_inatm200stn_processed_data_session01.csv'); ... 
                                          fullfile('Dataset_C_20181024','2018-10-24_Session02','02_processed','gnss_inatm200stn_processed_data_session02.csv'); ... 
                                          fullfile('Dataset_C_20181024','2018-10-24_Session03','02_processed','gnss_inatm200stn_processed_data_session03.csv'); ... 
                                        };
                         
imu_inatm200stn_processed_data_paths = { ... 
                                         fullfile('Dataset_C_20181024','2018-10-24_Session01','02_processed','imu_inatm200stn_processed_data_session01.csv'); ...
                                         fullfile('Dataset_C_20181024','2018-10-24_Session02','02_processed','imu_inatm200stn_processed_data_session02.csv'); ...
                                         fullfile('Dataset_C_20181024','2018-10-24_Session03','02_processed','imu_inatm200stn_processed_data_session03_part01.csv'); ...
                                         fullfile('Dataset_C_20181024','2018-10-24_Session03','02_processed','imu_inatm200stn_processed_data_session03_part02.csv'); ...
                                         fullfile('Dataset_C_20181024','2018-10-24_Session03','02_processed','imu_inatm200stn_processed_data_session03_part03.csv'); ... 
                                       };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNSS: Thales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gnss_thales_limited_processed_data_root_string = ... 
    [ ... 
      gnss_thales_processed_data_root_string, ... 
      '__', ... 
      num2str(start_time([12,13,15,16,18,19])), ...              
      '_', ... 
      num2str(end_time([12,13,15,16,18,19])) ... 
    ];

[out_data, num_files, load_flag] = importFromFile(gnss_thales_processed_data_paths,gnss_thales_limited_processed_data_root_string,'DataType','processed','StartTime',start_time,'EndTime',end_time);
for i = find( (load_flag(:)'==1) | (load_flag(:)'==3) ) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_thales_limited_processed_data_root_string,'_session',sprintf('%02i',i)],gnss_thales_limited_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_thales_processed_data_root_string,out_data); 
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNSS: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gnss_inatm200stn_limited_processed_data_root_string = ... 
    [ ... 
      gnss_inatm200stn_processed_data_root_string, ... 
      '__', ... 
      num2str(start_time([12,13,15,16,18,19])), ...              
      '_', ... 
      num2str(end_time([12,13,15,16,18,19])) ... 
    ];

[out_data, num_files, load_flag] = importFromFile(gnss_inatm200stn_processed_data_paths,gnss_inatm200stn_limited_processed_data_root_string,'DataType','processed','StartTime',start_time,'EndTime',end_time);
for i = find( (load_flag(:)'==1) | (load_flag(:)'==3) ) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[gnss_inatm200stn_limited_processed_data_root_string,'_session',sprintf('%02i',i)],gnss_inatm200stn_limited_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',gnss_inatm200stn_processed_data_root_string,out_data); 
clear out_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU: iNat M200 STN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imu_inatm200stn_limited_processed_data_root_string = ... 
    [ ... 
      imu_inatm200stn_processed_data_root_string, ... 
      '__', ... 
      num2str(start_time([12,13,15,16,18,19])), ...              
      '_', ... 
      num2str(end_time([12,13,15,16,18,19])) ... 
    ];

[out_data, num_files, load_flag] = importFromFile(imu_inatm200stn_processed_data_paths,imu_inatm200stn_limited_processed_data_root_string,'DataType','processed','StartTime',start_time,'EndTime',end_time);
for i = find( (load_flag(:)'==1) | (load_flag(:)'==3) ) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[imu_inatm200stn_limited_processed_data_root_string,'_session',sprintf('%02i',i)],imu_inatm200stn_limited_processed_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',imu_inatm200stn_processed_data_root_string,out_data);
clear out_data

%% Finish script

import_C_limited_processed_data_executed = true;
