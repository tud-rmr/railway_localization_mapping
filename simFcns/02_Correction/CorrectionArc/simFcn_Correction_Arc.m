function [s_hat_arc,P_s_hat_arc] = simFcn_Correction_Arc(s_0,P_s_0,z_s,R_s,Q,reset_flag)
% [p_hat_arc,P_p_hat_arc] = simFcn_Correction_Arc(s_0,P_s_0,z_s,R_s,Q,reset_flag)
%

coder.extrinsic('calcDiscreteKalmanFilter','isPositiveDefinite','nearestSPD');



% Initialization and checks _______________________________________________

persistent mem_var_s_0 mem_var_P_s_0 

if isempty(mem_var_s_0) || isempty(mem_var_P_s_0) || reset_flag
    mem_var_s_0 = s_0;    
    mem_var_P_s_0 = P_s_0; 
    
    s_hat_arc = mem_var_s_0;
    P_s_hat_arc = mem_var_P_s_0;

    return
end % if

% % System convariance matrix
% Q = zeros(length(mem_var_s_0)); 
% % Q = P_s_0;
% % Q = blkdiag(diag([20 20]),zeros(6));
% Q = blkdiag(zeros(2),0.01*0.1,0.01*0.1,deg2rad(0.001*0.1),0.01*0.1,0,0);

% Plug-in states that have been estimated outside
mem_var_s_0(7) = s_0(7);
mem_var_s_0(8) = s_0(8);

P_0 = blkdiag(mem_var_P_s_0(1:6,1:6),P_s_0(7:8,7:8));
if isPositiveDefinite(P_0,1e-9)
    mem_var_P_s_0 = P_0;
else
    mem_var_P_s_0 = nearestSPD(P_0);
end % if

% Kalman filtering ________________________________________________________
%   x = [x, y, x_0, y_0, psi_0, r, l, w]
%   z = [x, y, v, psi];

[x_temp_hat,P_x_temp_hat,~,~] = calcDiscreteKalmanFilter([],z_s(1:3),Q,R_s(1:3,1:3),mem_var_s_0,mem_var_P_s_0,'simFcn_Correction_Arc_StateTransitionFcn','simFcn_Correction_Arc_MeasurementFcn',[]);
% x_temp_hat(5) = deg2rad(simplifyHeadingD(rad2deg(x_temp_hat(5)),'180'));
mem_var_s_0 = x_temp_hat;
mem_var_P_s_0 = P_x_temp_hat;

% Output __________________________________________________________________

s_hat_arc = x_temp_hat;
P_s_hat_arc = P_x_temp_hat;

end % end function


