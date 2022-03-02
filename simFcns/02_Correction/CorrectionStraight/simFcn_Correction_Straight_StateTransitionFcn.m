function [x_hat,F] = simFcn_Correction_Straight_StateTransitionFcn(x_hat_old,z,u,sample_time)
%  [x_hat,F] = simFcn_Correction_Straight_StateTransitionFcn(x_hat_old,z,u,sample_time)

x     = x_hat_old(1);
y     = x_hat_old(2);
x_0   = x_hat_old(3);
y_0   = x_hat_old(4);
psi_0 = x_hat_old(5);
l     = x_hat_old(6);

x_hat(1,1) = l*sin(psi_0) + x_0;
x_hat(2,1) = l*cos(psi_0) + y_0;
x_hat(3,1) = x_0;
x_hat(4,1) = y_0;
x_hat(5,1) = psi_0;
x_hat(6,1) = l;

F = [ 0, 0, 1, 0,   l*cos(psi_0), sin(psi_0);... 
      0, 0, 0, 1,  -l*sin(psi_0), cos(psi_0); ... 
      0, 0, 1, 0,             0,          0; ... 
      0, 0, 0, 1,             0,          0; ... 
      0, 0, 0, 0,             1,          0; ... 
      0, 0, 0, 0,             0,          1];

end % end function 

