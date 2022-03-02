function [geometry,reset_flag,s_0,P_s_0,z_s,R_s,l_old,gamma] = simFcn_TgcDataPreparation(enable,x_imm,P_imm,mu,z,R,d_bogies,v_min_tgc)
% [geometry,reset_flag,s_0,P_s_0,z_s,R_s,l_old,gamma] = simFcn_TgcDataPreparation(enable,x_imm,P_imm,mu,z,R,d_bogies)
% 

%% Init

persistent mem_var_v_low mem_var_geometry mem_var_d mem_var_P_d

if isempty(mem_var_v_low)
    mem_var_v_low = 0;
end % if

if isempty(mem_var_geometry)
    mem_var_geometry = 0;
end % if

if isempty(mem_var_d)
    mem_var_d = 0;
end % if

if isempty(mem_var_P_d)
    mem_var_P_d = P_imm(3,3);
end % if

l_old = nan(2,1);

% Decompose IMM state vector
x   = x_imm(1);
y   = x_imm(2);
d   = x_imm(3);
v   = x_imm(4);
a   = x_imm(5);
psi = x_imm(6);
w   = x_imm(7);

%% Calculation: Geometry detection / change detection

[~,max_mu_index] = max(mu);
switch max_mu_index
    case 1 % straight _____________________________________________________
        geometry = 1;
    case 2 % arc __________________________________________________________
        geometry = 3;
    otherwise 
        geometry = 0;
end % switch

% Disabled from outside
if ~enable
    geometry = 0;
end % if

% Too slow for meaningful correction
% min_speed = round(2*1/3.6 * 10e2)/10e2;
% switch_hysteresis = round(0.4/3.6 * 10e2)/10e2;
min_speed = v_min_tgc(1);
switch_hysteresis = v_min_tgc(2);
if abs(x_imm(4)) < (min_speed-switch_hysteresis/2)
    mem_var_v_low = 1;
    geometry = 0;
elseif abs(x_imm(4)) > (min_speed+switch_hysteresis/2)
    mem_var_v_low = 0;
elseif mem_var_v_low
    geometry = 0;    
end % if

% Check for new geometry
if (geometry ~= mem_var_geometry)
    reset_flag = 1;
else 
    reset_flag = 0;
end % if
mem_var_geometry = geometry;

%% Calculation: Set initial parameters

l_old(1) = d-mem_var_d;
l_old(2) = P_imm(3,3)+mem_var_P_d;

if reset_flag % new geometry detected
    
    mem_var_d = d;
    mem_var_P_d = P_imm(3,3);
    
    l = 0;
    var_l = mem_var_P_d;
        
else
    
%     l = d-mem_var_d;
%     var_l = P_imm(3,3)+mem_var_P_d;
    l = l_old(1);
    var_l = l_old(2);
    
end % if

switch geometry
    
    case 1 % straight _____________________________________________________
        
        gamma = 0;
        
        s_0 = [ x; y; x; y; psi; 0; l; 0 ];
        P_s_0 = blkdiag(P_imm(1:2,1:2),P_imm(1:2,1:2),P_imm(6,6),0,var_l,0);
        
    otherwise
        
        % Check for small turn rates
        if abs(w) < 1e-7
            
            %error('simFcn_TgcDataPreparation: No ''straight'' detected, altough ''w'' is very small!')
            gamma = 0;
            s_0 = [ x; y; x; y; psi; 0; l; w ];
            P_s_0 = blkdiag(P_imm(1:2,1:2),P_imm(1:2,1:2),P_imm(6,6),0,var_l,P_imm(7,7)); 

        else
            
            % Correct misalignment between track and passenger-cabin
            gamma = calcBogieVehicleMisalignment(v,w,d_bogies);
            dir = sign(v);
            psi = psi - dir*gamma;

            s_0 = [ x; y; x; y; psi; v/w; l; w ];
            P_s_0 = blkdiag(P_imm(1:2,1:2),P_imm(1:2,1:2),P_imm(6,6),P_imm(4,4)/w^2,var_l,P_imm(7,7)); 
            
        end % if
        
end % switch

z_s = [x; y; v; psi];
R_s = P_imm([1:2,4,6],[1:2,4,6]);

end % function

%% Helper functions

function gamma = calcBogieVehicleMisalignment(v,w,d_bogies)

% Defines _________________________________________________________________

min_turn_rate = 1e-2; % minimal detectable turn rate in rad/s
min_speed = 1 / 3.6; % minimal speed in m/s
min_radius = 30; % minimal realistic curve radius in m

% Calculations ____________________________________________________________

gamma = 0; % standard output if it is not overwritten below

if (abs(w) > min_turn_rate) && (abs(v) > min_speed)
    radius = v/w;
    
    if abs(radius) > min_radius
        gamma = -asin( d_bogies/(2*radius) );
    end % if 
end % if

end % function
