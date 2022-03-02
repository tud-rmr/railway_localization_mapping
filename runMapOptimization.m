% 
% Run optimization to create a continuous railway map
%   
%   Author: Hanno Winter
%   Date: 28-Mar-2021; Last revision: 03-Jul-2021

if ~exist('init_localization_executed','var')
    error('runMapOptimization.m: Initaliziation script not executed!')
end % if
if ~exist('simout_optimization_data_selector','var')
    error('runMapOptimization.m: SimOut data not available!')
end % if

clear prepare_optim_out_executed % for convience: clear this variable to make sure that the script to process the optimization results can be executed
clear map_data_prepared % for convience: clear this variable to make sure that the script to process the optimization results can be executed

%% Settings

optim_track_ids = []; % 240:287 % 80:130 % 50:90 % Limit railway-map to this range of IDs

%% Init

% Load config
switch input_data_selector
    
    case {'C'}
        
        configCOptimization
                
    case {'BS'}
        
        configBsOptimization
        
    case {'NT'}
        
        configNtOptimization
                        
end % switch

% Prepare data from Simulink
prepareSimOutData 

% Prepare data for optimization
prepareOptimInData 

%% Optimization
   
tic

z_start_test = [];
z_end_test = [];
for optim_i = 1:size(optim_start_indices,1)
    
    fprintf('\nOptimization Run: %i/%i\n',optim_i,size(optim_start_indices,1));

    % Init ________________________________________________________________
    
    % Select current optimization frame 
    map_indices = (optim_start_indices(optim_i):optim_end_indices(optim_i));                
    map_i.topology = optim_map.topology(map_indices,map_indices);
    map_i.track_start_points = optim_map.track_start_points(map_indices,:);
    map_i.track_maps = optim_map.track_maps(map_indices,:);
    optim_data_selector_i = optim_data_selector(map_indices,:);
    optim_data_selector_i(:,1) = (1:length(map_indices))';
    
    % Convert railway-map to "vector" representation
    [optim_x_0,optim_x_0_descriptor] = railwayMap2X(map_i);

    % Boundaries
    switch options.Algorithm
        case {'levenberg-marquardt'}
            lb = [];
            ub = [];
        otherwise
            lb = -inf(size(optim_x_0));
            for i = find(optim_x_0_descriptor(1,:) == 3)
                lb(1,i) = 0;
                lb(3,i) = 0;
                lb(4,i) = 0;
            end % for i    
            ub = inf(size(optim_x_0));
    end % switch

    % Optimization ________________________________________________________

    fun_handle = @(x) railwayMapErrorFcn( ... 
                                          x, ... 
                                          optim_x_0_descriptor, ... 
                                          optim_input_data, ... 
                                          optim_input_data_P, ... 
                                          optim_data_selector_i, ... 
                                          error_fcn_config ... 
                                        );
    [optim_x,~,optim_residual,~,optim_output] = ... 
        lsqnonlin(fun_handle,optim_x_0(:),lb(:),ub(:),options);
    
    % Processing __________________________________________________________
    
    % Reshape
    optim_x = reshape(optim_x,4,length(optim_x(:))/4);
    
    % Correct negative length
    optim_x = ensurePositiveLength(optim_x,optim_x_0_descriptor);
    
    % Correct huge radii
    % optim_x = ensureRadiusLimits(optim_x,optim_x_0_descriptor);

    % Convert to continous railway-map
    optim_map_temp_i = x2contTableRailwayMap(optim_x,optim_x_0_descriptor,[]);

    % Limit start and end of railway-map
    z_start = optim_input_data(:,logical(optim_data_selector(map_indices(1),2:end)));    
    if isempty(z_start)
        z_start = zeros(2,1);
    else
        z_start = z_start(:,1);
    end % if
    
    z_end = optim_input_data(:,logical(optim_data_selector(map_indices(end),2:end)));
    if ~isempty(z_end)            
        z_end = z_end(:,end);          
        if optim_i == 1 % || optim_i == size(optim_start_indices,1)
            optim_map_temp_i = limitRailwayMap(z_start,z_end,optim_map_temp_i);
        else
            optim_map_temp_i = limitRailwayMap([],z_end,optim_map_temp_i);
        end % if
    end % if

    % Update overall railway-map
    package_indices = (package_start_indices(optim_i):package_end_indices(optim_i));
    temp_i_indices = (package_start_indices(optim_i):package_end_indices(optim_i)) - optim_start_indices(optim_i)+1;
        
    optim_map.track_start_points(package_indices,:) = optim_map_temp_i.track_start_points(temp_i_indices,:);
    optim_map.track_start_points.ID = (1:size(optim_map.track_start_points,1))';
    optim_map.track_maps(package_indices,:) = optim_map_temp_i.track_maps(temp_i_indices,:);
    optim_map.track_maps.ID = (1:size(optim_map.track_maps,1))';
        
end % for

% Correct final optimized map _____________________________________________

optim_map.topology = zeros(size(optim_map.topology));
optim_map.track_maps.track_element((optim_map.track_maps.track_element == 11),:) = 1;

%% Finish

optim_time = toc;

run_optimization_executed = true;

fprintf('\nFinished optimization (elapsed time: %.0fs)\n\n',optim_time)
