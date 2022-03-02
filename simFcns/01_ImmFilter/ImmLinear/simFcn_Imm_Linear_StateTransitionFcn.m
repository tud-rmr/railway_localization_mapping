function [x_hat,F] = simFcn_Imm_Linear_StateTransitionFcn(x_hat_old,z,u,T)

% System __________________________________________________________________

% u1 = u(1);
% u2 = u(2);

x_hat = x_hat_old;

x   = x_hat_old(1);
y   = x_hat_old(2);
d   = x_hat_old(3);
v   = x_hat_old(4);
a   = x_hat_old(5);
% a   = u1;
psi = x_hat_old(6);
w   = x_hat_old(7);
% w   = u2;

% delta_x = 1/2*a*T^2*sin(psi) + v*T*sin(psi);
% delta_y = 1/2*a*T^2*cos(psi) + v*T*cos(psi);

delta_x = 1/2*a*T^2*sin(psi) + v*T*sin(psi);
delta_y = 1/2*a*T^2*cos(psi) + v*T*cos(psi);

x_hat(1) = x + delta_x;
x_hat(2) = y + delta_y;
x_hat(3) = d + v*T + 1/2*a*T^2;
x_hat(4) = v + a*T;
x_hat(5) = a;
x_hat(6) = psi;
x_hat(7) = 0;

% Linearization ___________________________________________________________

x   = x_hat(1);
y   = x_hat(2);
d   = x_hat(3);
v   = x_hat(4);
a   = x_hat(5);
psi = x_hat(6);
w   = x_hat(7);

F = [ 1, 0, 0, T*sin(psi), (T^2*sin(psi))/2,  (T*cos(psi)*(2*v + T*a))/2, 0; ... 
      0, 1, 0, T*cos(psi), (T^2*cos(psi))/2, -(T*sin(psi)*(2*v + T*a))/2, 0; ... 
      0, 0, 1,          T,            T^2/2,                           0, 0; ... 
      0, 0, 0,          1,                T,                           0, 0; ... 
      0, 0, 0,          0,                1,                           0, 0; ... 
      0, 0, 0,          0,                0,                           1, 0; ... 
      0, 0, 0,          0,                0,                           0, 0];

end % function

