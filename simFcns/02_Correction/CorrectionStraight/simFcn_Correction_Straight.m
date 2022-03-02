function [s_hat_straight,P_s_hat_straight] = simFcn_Correction_Straight(s_0,P_s_0,z_s,R_s,Q,reset_flag)
% [s_hat_straight,P_s_hat_straight] = simFcn_Correction_Straight(s_0,P_s_0,z_s,R_s,Q,reset_flag)
%

coder.extrinsic('calcDiscreteKalmanFilter','isPositiveDefinite','nearestSPD');

% Initialization and checks _______________________________________________

persistent mem_var_s_0 mem_var_P_s_0 

% Map general state vector to straight specific vector
x_indices = [1,2,3,4,5,7];
s_0 = s_0(x_indices);
P_s_0 = P_s_0(x_indices,x_indices);

% Map general measurement vector to straight specific vector
z_indices = [1,2,4];
z_s = z_s(z_indices);
R_s = R_s(z_indices,z_indices);

% (Re)-Initialize filter
if isempty(mem_var_s_0) || isempty(mem_var_P_s_0) || reset_flag
    mem_var_s_0 = s_0;    
    mem_var_P_s_0 = P_s_0;
    
    s_hat_straight = mem_var_s_0;
    P_s_hat_straight = mem_var_P_s_0;

    return
end % if

% % System convariance matrix
% Q = zeros(length(mem_var_s_0)); 
% % Q = P_s_0;
% % Q = blkdiag(diag([20 20]),zeros(4));
% Q = blkdiag(zeros(2),0.01*0.1,0.01*0.1,deg2rad(0.001*0.1),0);

% Plug-in states that have been estimated outside
mem_var_s_0(6) = s_0(6);

P_0 = blkdiag(mem_var_P_s_0(1:5,1:5),P_s_0(6,6));
% P_0 = mem_var_P_s_0;
% P_0(6,6) = P_s_0(6,6);
if isPositiveDefinite(P_0,1e-9)
    mem_var_P_s_0 = P_0;
else
    mem_var_P_s_0 = nearestSPD(P_0);
end % if

% Kalman filtering ________________________________________________________
%   x = [x, y, x_0, y_0, psi_0, l]
%   z = [x, y, psi];

[x_temp_hat,P_x_temp_hat,~,~] = calcDiscreteKalmanFilter([],z_s(1:2),Q,R_s(1:2,1:2),mem_var_s_0,mem_var_P_s_0,'simFcn_Correction_Straight_StateTransitionFcn','simFcn_Correction_Straight_MeasurementFcn',[]);
% x_temp_hat(5) = deg2rad(simplifyHeadingD(rad2deg(x_temp_hat(5)),'180'));
mem_var_s_0 = x_temp_hat;
mem_var_P_s_0 = P_x_temp_hat;

% Output __________________________________________________________________

s_hat_straight = x_temp_hat;
P_s_hat_straight = P_x_temp_hat;

end % end function
