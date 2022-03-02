% 
% Prepare output data from map creation via optimization
%   
%   Author: Hanno Winter
%   Date: 08-Apr-2021; Last revision: 24-Apr-2021

%% Init

if ~exist('prepare_optim_out_executed','var') || ~prepare_optim_out_executed
    prepare_optim_out_executed = false;
    fprintf('Processing optimization output data\n');
elseif prepare_optim_out_executed
    fprintf('Already processed optimization output data\n');
    return  
end % if

prepareSimOutData

%% Calculate Positions

% Calculate railway-map positions
optim_map_density = 1; % density in m
[~,~,optim_map_abs_pos,t_optim_map_abs_pos,~,~,~,~] = calcMapProperties(optim_map,optim_map_density);

% Convert railway-map postions to lat-lon coordinates
[optim_map_abs_position_latitude,optim_map_abs_position_longitude] = ... 
    utm2ll(p_0_utm(1)+optim_map_abs_pos(1,:),p_0_utm(2)+optim_map_abs_pos(2,:),p_0_utm(3),'wgs84');


%% Finish

prepare_optim_out_executed = true;

%% Save

prompt_str = 'Do you want to save the processed optimout data (j/n)?';
user_str = input(prompt_str,'s');

switch user_str
    case {'j','y','yes','ja'}
        optimout_data_saved = false;
        save_optim_out_data = true;
    case {'n','no','nein'}
        save_optim_out_data = false;
    otherwise
        save_optim_out_data = false;
end % switch

if save_optim_out_data
    saveOptimData
end %if
