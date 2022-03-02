function [x_hat,F] = simFcn_Imm_Trash_StateTransitionFcn(x_hat_old,z,u,T)

% [x_hat,F] = simFcn_Imm_Circle_StateTransitionFcn(x_hat_old,z,u,T);

% Switch to CA model for small turn rates _________________________________

% if abs(x_hat_old(7)) > 5*eps
% if (abs(u(2)) < 1e-6) || (abs(x_hat_old(7)) < 1e-9)
% if (abs(z(4)) > 0.1*(pi/180))
if (abs(x_hat_old(7)) < 1e-7)

    [x_hat,F] = simFcn_Imm_Linear_StateTransitionFcn(x_hat_old,z,u,T);
    return
    
end % if

% if (abs(x_hat_old(7)) < 0.005) && sign(x_hat_old(7)) ~= 0
% 
%     x_hat_old(7) = sign(x_hat_old(7)) * 0.005;
%     
% elseif (abs(x_hat_old(7)) < 0.005) && sign(x_hat_old(7)) == 0
%     
%     x_hat_old(7) = 0.005;
%     
% end % if

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

% delta_x = 1/w^2*((-v*w - a*w*T)*cos(psi+w*T) ... 
%           + a*sin(psi+w*T)                   ... 
%           + v*w*cos(psi) - a*sin(psi));
% delta_y = 1/w^2*((v*w + a*w*T)*sin(psi+w*T)  ... 
%           + a*cos(psi+w*T)                   ... 
%           - v*w*sin(psi) - a*cos(psi));  
      
delta_x = 1/w^2*((-v*w - a*w*T)*cos(psi+w*T) ... 
          + a*sin(psi+w*T)                   ... 
          + v*w*cos(psi) - a*sin(psi));
delta_y = 1/w^2*((v*w + a*w*T)*sin(psi+w*T)  ... 
          + a*cos(psi+w*T)                   ... 
          - v*w*sin(psi) - a*cos(psi));  
      
x_hat(1,1) = x + delta_x;
x_hat(2,1) = y + delta_y;
x_hat(3,1) = d + v*T + 1/2*a*T^2;
x_hat(4,1) = v + a*T;
x_hat(5,1) = a;
x_hat(6,1) = psi + w*T;
x_hat(7,1) = w;

% Linearization ___________________________________________________________

x   = x_hat(1);
y   = x_hat(2);
d   = x_hat(3);
v   = x_hat(4);
a   = x_hat(5);
psi = x_hat(6);
w   = x_hat(7);

F = [ 1, 0, 0, -(cos(psi + T*w) - cos(psi))/w, -(sin(psi) - sin(psi + T*w) + T*w*cos(psi + T*w))/w^2,  (a*cos(psi + T*w) - a*cos(psi) + w*sin(psi + T*w)*(v + T*a) - v*w*sin(psi))/w^2, (2*a*sin(psi) - 2*a*sin(psi + T*w) - v*w*cos(psi) + v*w*cos(psi + T*w) + T*v*w^2*sin(psi + T*w) + T^2*a*w^2*sin(psi + T*w) + 2*T*a*w*cos(psi + T*w))/w^3; ...
      0, 1, 0,  (sin(psi + T*w) - sin(psi))/w,  (cos(psi + T*w) - cos(psi) + T*w*sin(psi + T*w))/w^2, -(a*sin(psi + T*w) - a*sin(psi) - w*cos(psi + T*w)*(v + T*a) + v*w*cos(psi))/w^2, (2*a*cos(psi) - 2*a*cos(psi + T*w) + v*w*sin(psi) - v*w*sin(psi + T*w) + T*v*w^2*cos(psi + T*w) + T^2*a*w^2*cos(psi + T*w) - 2*T*a*w*sin(psi + T*w))/w^3; ...
      0, 0, 1,                              T,                                                 T^2/2,                                                                                0,                                                                                                                                                        0; ...
      0, 0, 0,                              1,                                                     T,                                                                                0,                                                                                                                                                        0; ...
      0, 0, 0,                              0,                                                     1,                                                                                0,                                                                                                                                                        0; ...
      0, 0, 0,                              0,                                                     0,                                                                                1,                                                                                                                                                        T; ...
      0, 0, 0,                              0,                                                     0,                                                                                0,                                                                                                                                                        1];
 
  
end % function
