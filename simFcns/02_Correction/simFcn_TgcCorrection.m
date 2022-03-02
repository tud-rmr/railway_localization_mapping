function [s,P_s,o,P_o] = simFcn_TgcCorrection(enable,reset_flag,x,P_x,s_0,P_s_0,z_s,R_s,Q_straight,Q_arc,geometry)
% [s,P_s,o,P_o] = simFcn_TgcCorrection(enable,reset_flag,x,P_x,s_0,P_s_0,z_s,R_s,geometry)
%

%% Calculations

if ~enable % this condition is acutally not necessary, because it should have been checked before.
    geometry = 0;
end % if

switch geometry

    case 1 % straight _____________________________________________________

        [s,P_s] = simFcn_Correction_Straight(s_0,P_s_0,z_s,R_s,Q_straight,reset_flag);
        s = [s(1:5);0;s(6);0];
        P_s = blkdiag(P_s(1:5,1:5),0,P_s(6,6),0);

    case 3 % arc __________________________________________________________

        [s,P_s] = simFcn_Correction_Arc(s_0,P_s_0,z_s,R_s,Q_arc,reset_flag);

    otherwise % ___________________________________________________________
        
        %s = s_0;
        %P_s = P_s_0;        
        s = [s_0(1:5);0;s_0(7:8)];
        P_s = blkdiag(P_s_0(1:5,1:5),0,P_s_0(7:8,7:8));

end % switch

o = [s(1:2);x(4);x(6)];
P_o = blkdiag(P_s(1:2,1:2),P_x(4,4),P_x(6,6));

end % function
