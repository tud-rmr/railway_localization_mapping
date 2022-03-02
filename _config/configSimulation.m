% 
% Configuration file for simulation
%   
%   Author: Hanno Winter
%   Date: 18-Jul-2021; Last revision: 18-Jul-2021

%%

% Dataset _________________________________________________________________
input_data_selector = 'BS'; % Choose dataset {'C','BS','NT'}

% Simulink model __________________________________________________________
sim_model_name = 'localization';

% Artificial GPS outages __________________________________________________
insert_gps_outages = false; % enable or disable additional integration of GPS outages
target_gps_availability = 0.75; % target GPS availability [0,1]
target_gps_outage_length = [1 20]; % length of additional GPS outages as limits of a uniform distribution
 