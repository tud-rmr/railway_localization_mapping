function [r_circle,r_length,omega_r_circle,p_circle_ppd,l_circle_ppd] = calcErrorToCircularGeometry(z,R,z_data_selector,track_maps,track_start_points)
% [r_circle,omega_r_circle,p_circle_ppd,l_circle_ppd] = calcErrorToCircularGeometry(z,R,z_data_selector,track_maps,track_start_points)
% 

coder.extrinsic('invChol_mex')

%% Init

tm_id_idx = 1;
tm_element_idx = 2;
tm_rend_idx = 4;
tm_length_idx = 5;
sp_id_idx = 1;
sp_x0_idx = 2;
sp_y0_idx = 3;
sp_phi0_idx = 4;

circle_selector = (track_maps(:,tm_element_idx) == 3);
circle_track_ids = track_maps(circle_selector,tm_id_idx);
circle_track_indices = find(ismember(track_start_points(:,sp_id_idx),circle_track_ids));

%% Calc

n = length(circle_track_ids);
r_circle = cell(n,1);
r_length = cell(n,1);
p_circle_ppd = cell(n,1);
l_circle_ppd = cell(n,1);
omega_r_circle = cell(n,1);
for i = 1:n
    current_track_index = circle_track_indices(i);
    current_track_ID = track_start_points(current_track_index,sp_id_idx);
    
    % Measurement data selection __________________________________________
    z_i_selector = logical(z_data_selector(z_data_selector(:,1) == current_track_ID,2:end));
    z_i = z(:,z_i_selector);
    R_i = R(:,:,z_i_selector);
    
    if isempty(z_i) % || size(z_i,2) < 2
        continue;
    end % if    
    
    % Calculation of fixed auxilary variables _____________________________
    p_0 = [ track_start_points(current_track_index,sp_x0_idx); ... 
            track_start_points(current_track_index,sp_y0_idx)];
    phi_0 = track_start_points(current_track_index,sp_phi0_idx);
    radius = track_maps(current_track_index,tm_rend_idx) * -1; % convert from railway convention to mathematical convention                   
    
    if radius >= 0        
        p_center = p_0 + abs(radius)*[-sind(phi_0);cosd(phi_0)]; 
    else
        p_center = p_0 + abs(radius)*[sind(phi_0);-cosd(phi_0)];              
    end % if
    
    z_hat_i = NaN(2,size(z_i,2));
    omega_r_arc_i = NaN(2,2,size(z_i,2));
    
    p_circle_ppd_i = nan(2,size(z_i,2));
    l_circle_ppd_i = nan(1,size(z_i,2));
    for j = 1:size(z_i,2)
        % z_hat ___________________________________________________________
        p_test = z_i(:,j);        
        p_00 = p_center + abs(radius)*(p_test-p_center)/norm((p_test-p_center));
%         p_circle_ppd_i = [p_circle_ppd_i,p_00];
        p_circle_ppd_i(:,j) = p_00;
        z_hat_i(:,j) = p_00;
        
        % Calulation of arc length ________________________________________
        a_vec = p_0-p_center;  
        b_vec = p_00-p_center;
        phi_a = mod(atan2d(a_vec(2),a_vec(1)),360);
        phi_b = mod(atan2d(b_vec(2),b_vec(1)),360);
        
        % reference all angles to vector a
        if radius > 0
            phi_b = mod(phi_b-phi_a,360);
            phi_a = 0;
        end % if
        if radius < 0
            phi_b = mod(phi_a-phi_b,360);
            phi_a = 0;
        end % if
        gamma = phi_b;
        if gamma > 180
            gamma = -(360-gamma);
        end % if        
%         l_circle_ppd_i = [l_circle_ppd_i,gamma/360*2*pi*abs(radius)];
        l_circle_ppd_i(:,j) = gamma/360*2*pi*abs(radius);
        
        
        % omega ___________________________________________________________
        % omega_r_arc_i(:,:,j) = R_i(:,:,j)^-1;
        if det(R_i(:,:,j)) == 0
            omega_r_arc_i(:,:,j) = inf(2);
        else
            omega_r_arc_i(:,:,j) = invChol_mex(R_i(:,:,j));
        end % if 
    end % for j
    p_circle_ppd{i} = p_circle_ppd_i;
    l_circle_ppd{i} = l_circle_ppd_i;
    
%         % TEST %%%
%     plot(z_i(1,:),z_i(2,:),'bx','MarkerSize',10); hold on;
% %     L = track_maps(current_track_index,tm_length_idx);
% %     l = linspace(0,L,30);
% %     p_arc_x = radius*(1-cosd(rad2deg(l./radius)+phi_0))+p_0(1)-radius*(1-cosd(phi_0));
% %     p_arc_y = radius*(  sind(rad2deg(l./radius)+phi_0))+p_0(2)-radius*(  sind(phi_0));
% %     p_arc = [p_arc_x;p_arc_y];
% %     plot(p_arc(1,:),p_arc(2,:),'-','LineWidth',1.5);
%     plot(p_circle_ppd{i}(1,:),p_circle_ppd{i}(2,:),'-','LineWidth',1.5);
%     axis equal
%     % ^^^ %%%
    
    % Ouput for track 'i' _________________________________________________
    r_circle{i} = z_i - z_hat_i;
    r_length{i} = [ norm(z_i(:,1)-p_circle_ppd{i}(:,1)); norm(z_i(:,end)-p_circle_ppd{i}(:,end)) ];
    omega_r_circle{i} = omega_r_arc_i;
end % for i

end % function
