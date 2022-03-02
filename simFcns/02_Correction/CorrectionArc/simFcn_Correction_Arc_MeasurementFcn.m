function [z_hat,H] = simFcn_Correction_Arc_MeasurementFcn(x_hat,z,u,sample_time)
% [z_hat,H] = simFcn_Correction_Arc_MeasurementFcn(x_hat,u,sample_time)
%

x     = x_hat(1);
y     = x_hat(2);
x_0   = x_hat(3);
y_0   = x_hat(4);
psi_0 = x_hat(5);
r     = x_hat(6);
l     = x_hat(7);
w     = x_hat(8);

z_hat(1,1) = x;
z_hat(2,1) = y;
z_hat(3,1) = w*r;
% z_hat(4,1) = l/r+psi_0;

% H = [ 1, 0, 0, 0, 0,      0,   0, 0;   
%       0, 1, 0, 0, 0,      0,   0, 0;   
%       0, 0, 0, 0, 0,      w,   0, r;   
%       0, 0, 0, 0, 1, -l/r^2, 1/r, 0];
  
H = [ 1, 0, 0, 0, 0,      0,   0, 0;   
      0, 1, 0, 0, 0,      0,   0, 0;   
      0, 0, 0, 0, 0,      w,   0, r];

end % end function

