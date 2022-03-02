% 
% Initialize Track-Geometry-Constraint Correction of IMM data
%   
%   Author: Hanno Winter
%   Date: 25-Mar-2021; Last revision: 17-Apr-2021

%% Init

if ~exist('init_localization_executed','var')
    warning('initTgcCorrection.m: This script is supposed to be called from ''initLocalization.m''!')
    % return
% elseif ~exist('prepare_sim_input_data_executed','var')
%     warning('initMapBuilding.m: Run ''prepareSimInputData.m'' before executing this script!')
%     return
else
    init_tgc_correction_executed = false;
end % if

%% Calculations

Q_straight = blkdiag(0,0,sigma_p_0.^2,sigma_p_0.^2,sigma_psi_straight.^2,0);
Q_arc = blkdiag(0,0,sigma_p_0.^2,sigma_p_0.^2,sigma_psi_arc.^2,sigma_r.^2,0,0);

%% Finish

init_tgc_correction_executed = true;
