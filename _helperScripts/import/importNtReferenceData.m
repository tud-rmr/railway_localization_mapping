%
%   Data available at: http://dx.doi.org/10.25534/tudatalib-359
%
%   Author: Hanno Winter
%   Date: 09-Dec-2020; Last revision: 18-Apr-2021
% 

%% Init

loadNtFilenamePatterns

%% Checks

% if ~exist('import_NT_reference_data_executed','var')
%     import_NT_reference_data_executed = false;
% elseif import_NT_reference_data_executed
%     fprintf('Import reference data script already executed!\n');
%     return  
% end % if

%% Processed Data
%
%   Filepath or filename has to contain session information, i.e. it has to 
%   contain a string which can be found with the regular 
%   expression '[Ss]ession\d{2}'. Valid strings are, 
%   e.g. 'Session01' or 'session01'.
%
%   If the data has been chunked into several parts it is import to provide
%   the order of the parts.
%   

% ref_inatm200stn_data_paths = { ... 
%                                fullfile('ref_inatm200stn_internal_ekf_data_session01.csv'); ... 
%                              };
                         
ref_ekfFusion_inatm200stn_data_paths = { ... 
                                         fullfile('Dataset_NT_20190311','2019-03-11_Session01','03_reference','ref_ekfFusion_inatm200stn_data_session01.csv'); ... 
                                         fullfile('Dataset_NT_20190311','2019-03-11_Session02','03_reference','ref_ekfFusion_inatm200stn_data_session02.csv'); ...
                                         fullfile('Dataset_NT_20190311','2019-03-15_Session03','03_reference','ref_ekfFusion_inatm200stn_data_session03.csv'); ...
                                         fullfile('Dataset_NT_20190311','2019-03-15_Session04','03_reference','ref_ekfFusion_inatm200stn_data_session04.csv'); ...
                                         fullfile('Dataset_NT_20190311','2019-03-15_Session05','03_reference','ref_ekfFusion_inatm200stn_data_session05.csv') ...
                                       };
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference: EKF Fusion (iNat M200 STN GNSS and IMU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out_data, num_files, load_flag] = importFromFile(ref_ekfFusion_inatm200stn_data_paths,ref_ekfFusion_inatm200stn_reference_data_root_string,'DataType','reference');
for i = find(load_flag(:)'==1) % save to .mat for faster access in the future
    exportToFile(out_data(i,:),[ref_ekfFusion_inatm200stn_reference_data_root_string,'_session',sprintf('%02i',i)],ref_ekfFusion_inatm200stn_reference_data_root_string,'SaveTo','mat','NumFiles',num_files(i));
end % for i
assignin('base',ref_ekfFusion_inatm200stn_reference_data_root_string,out_data); 
clear out_data


%% Finish script

import_NT_reference_data_executed = true;
