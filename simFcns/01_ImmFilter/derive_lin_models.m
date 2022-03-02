clear all
close all
clc

syms x y d v a psi w T real
x_vec = [x,y,d,v,a,psi,w]';

%% Straight

x_hat_str(1) = x + 1/2*a*T^2*sin(psi) + v*T*sin(psi);
x_hat_str(2) = y + 1/2*a*T^2*cos(psi) + v*T*cos(psi);
x_hat_str(3) = d + v*T + 1/2*a*T^2;
x_hat_str(4) = v + a*T;
x_hat_str(5) = a;
x_hat_str(6) = psi;
x_hat_str(7) = 0;

F_straight = simplify(jacobian(x_hat_str,x_vec))

%% Circular arc

delta_x = 1/w^2*((-v*w - a*w*T)*cos(psi+w*T) ... 
          + a*sin(psi+w*T)                   ... 
          + v*w*cos(psi) - a*sin(psi));
delta_y = 1/w^2*((v*w + a*w*T)*sin(psi+w*T)  ... 
          + a*cos(psi+w*T)                   ... 
          - v*w*sin(psi) - a*cos(psi));   
      
x_hat_circ(1,1) = x + delta_x;
x_hat_circ(2,1) = y + delta_y;
x_hat_circ(3,1) = d + v*T + 1/2*a*T^2;
x_hat_circ(4,1) = v + a*T;
x_hat_circ(5,1) = a;
x_hat_circ(6,1) = psi + w*T;
x_hat_circ(7,1) = w;

F_circ = simplify(jacobian(x_hat_circ,x_vec))


