function  [x_hat,F] = simFcn_Correction_Arc_StateTransitionFcn(x_hat_old,z,u,sample_time)
% [x_hat,F] = simFcn_Correction_Arc_StateTransitionFcn(x_hat_old,u,sample_time)
%

x     = x_hat_old(1);
y     = x_hat_old(2);
x_0   = x_hat_old(3);
y_0   = x_hat_old(4);
psi_0 = x_hat_old(5);
r     = x_hat_old(6);
l     = x_hat_old(7);
w     = x_hat_old(8);

x_hat(1,1) = r*(1-cos(l/r+psi_0)) + (x_0 - r*(1-cos(psi_0)));
x_hat(2,1) = r*(  sin(l/r+psi_0)) + (y_0 - r*(  sin(psi_0)));
x_hat(3,1) = x_0;
x_hat(4,1) = y_0;
x_hat(5,1) = psi_0;
x_hat(6,1) = r;
x_hat(7,1) = l;
x_hat(8,1) = w;

F = [ 0, 0, 1, 0, r*(sin(psi_0 + l/r) - sin(psi_0)), cos(psi_0) - cos(psi_0 + l/r) - (l*sin(psi_0 + l/r))/r, sin(psi_0 + l/r), 0; ... 
      0, 0, 0, 1, r*(cos(psi_0 + l/r) - cos(psi_0)), sin(psi_0 + l/r) - sin(psi_0) - (l*cos(psi_0 + l/r))/r, cos(psi_0 + l/r), 0; ... 
      0, 0, 1, 0,                                 0,                                                      0,                0, 0; ... 
      0, 0, 0, 1,                                 0,                                                      0,                0, 0; ... 
      0, 0, 0, 0,                                 1,                                                      0,                0, 0; ... 
      0, 0, 0, 0,                                 0,                                                      1,                0, 0; ... 
      0, 0, 0, 0,                                 0,                                                      0,                1, 0; ... 
      0, 0, 0, 0,                                 0,                                                      0,                0, 1];

end % end function
