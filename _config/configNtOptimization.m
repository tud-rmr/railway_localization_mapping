% 
% Configuration file for map creation via optimization with the data set 
% 'Nürtingen' (NT) (see: https://tudatalib.ulb.tu-darmstadt.de/handle/tudatalib/2529)
%   
%   Author: Hanno Winter
%   Date: 08-Apr-2021; Last revision: 09-Oct-2021

%% Input data settings

optim_input_selector = 'TGC'; % Choose from {'GPS','IMM','TGC'}
z_data_dividor = 1; % Reduce measurement by data dividor
frame_size = 4; % Optimization frame size 
min_z = 20; % Minimal number of measurements for straights limiting the optimization frames
min_z_per_element.straight = 2; % Dismiss straights during optimization with less measurements
min_z_per_element.arc = 2; % Dismiss circular-arcs during optimization with less measurements
straight_sloppiness.threshold = 40; % Start punishing length of straights when it is getting longer than this value (in m), compared to the initial length
straight_sloppiness.range = 20; % Ramp length around straight sloppiness threshold to smooth punishing too long straights

%% Optimizer settings

options = optimoptions('lsqnonlin','Display','iter-detailed');
options.FiniteDifferenceStepSize = 1e-4;
options.ScaleProblem = 'jacobian';    
% options.FunctionTolerance = 1e-4;
options.StepTolerance = 1e-5;    
% options.MaxFunctionEvaluations = 100; 
% options.PlotFcn = 'myOptimPlotFcn'; % optimplotfval,optimplotstepsize,optimplotresnorm,myOptimPlotFcn
options.Algorithm = 'levenberg-marquardt';
options.UseParallel = true;
options.MaxIterations = 60;

error_fcn_config.min_z_per_element = min_z_per_element;
error_fcn_config.straight_sloppiness = straight_sloppiness;
