% 
% Initialize Localization
%   
%   Author: Hanno Winter
%   Date: 13-Mar-2021; Last revision: 19-Jul-2021

clear functions % Important: clear persistent variables!

%% Init

if ~exist('init_localization_executed','var')
    init_localization_executed = false;
    fprintf('Init localization\n');
elseif init_localization_executed
    fprintf('Already initialized localization\n');
    return  
end % if

%% Simulation Settings

% Load simulation config
configSimulation

%% Load and prepare input data

switch input_data_selector
    case {'C'}
        
        fprintf('\nUsing dataset: C -- Chemnitz\n\n')
        
        % Load config
        configCLocalization
        
        % Load data
        initCData
        
        % Prepare data
        prepareSimInData
        
        % Init Modules
        initStrapdown
        initImmFilter
        initTgcCorrection
        
    case {'BS'}
        
        fprintf('\nUsing dataset: BS -- Braunschweig\n\n')
        
        % Load config
        configBsLocalization
        
        % Load data
        initBsData
        
        % Prepare data
        prepareSimInData
        
        % Init Modules
        initStrapdown
        initImmFilter
        initTgcCorrection
        
    case {'NT'}
        
        fprintf('\nUsing dataset: NT -- Nürtingen\n\n')
        
        % Load config
        configNtLocalization
        
        % Load data
        initNtData
        
        % Prepare data
        prepareSimInData
        
        % Init Modules
        initStrapdown
        initImmFilter
        initTgcCorrection
        
    otherwise
        
        error('''initLocalization.m'' unsopported case!')
        
end % switch

%% Finish script

init_localization_executed = true;
