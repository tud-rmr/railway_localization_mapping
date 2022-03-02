function [z_hat,H] = simFcn_Correction_Straight_MeasurementFcn(x_hat,z,u,sample_time)
% [z_hat,H] = simFcn_Correction_Straight_MeasurementFcn(x_hat,z,u,sample_time)
%

x     = x_hat(1);
y     = x_hat(2);
x_0   = x_hat(3);
y_0   = x_hat(4);
psi_0 = x_hat(5);
l     = x_hat(6);

z_hat(1,1) = x;
z_hat(2,1) = y;
% z_hat(3,1) = psi_0;

% H = [ 1, 0, 0, 0, 0, 0; ... 
%       0, 1, 0, 0, 0, 0; ... 
%       0, 0, 0, 0, 1, 0];
  
H = [ 1, 0, 0, 0, 0, 0; ... 
      0, 1, 0, 0, 0, 0];

end % end function

