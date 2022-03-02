function [r_straight,r_length,omega_r_straight,p_straight_ppd,l_straight_ppd] = calcErrorToStraightGeometry(z,R,z_data_selector,track_maps,track_start_points)
% [r_straight,r_length,omega_r_straight,p_straight_ppd,l_straight_ppd] = calcErrorToStraightGeometry(z,R,z_data_selector,track_maps,track_start_points)
%

coder.extrinsic('invChol_mex')

%% Init

tm_id_idx = 1;
tm_element_idx = 2;
tm_length_idx = 5;
sp_id_idx = 1;
sp_x0_idx = 2;
sp_y0_idx = 3;
sp_phi0_idx = 4;

straight_selector = (track_maps(:,tm_element_idx) == 1);
% straight_selector = (track_maps(:,tm_element_idx) == 1) | (track_maps(:,tm_element_idx) == 11);
straight_track_ids = track_maps(straight_selector,tm_id_idx);
straigt_track_indices = find(ismember(track_start_points(:,sp_id_idx),straight_track_ids));

%% Calc

n = length(straight_track_ids);
r_straight = cell(n,1);
r_length = cell(n,1);
p_straight_ppd = cell(n,1);
l_straight_ppd = cell(n,1);
omega_r_straight = cell(n,1);
for i = 1:n
    current_track_index = straigt_track_indices(i);
    current_track_ID = track_start_points(current_track_index,sp_id_idx);
    
    p_0 = [track_start_points(current_track_index,sp_x0_idx);track_start_points(current_track_index,sp_y0_idx)];
    phi_0 = track_start_points(current_track_index,sp_phi0_idx);
    L = track_maps(current_track_index,tm_length_idx);
    t = [cosd(phi_0);sind(phi_0)];
    
    z_i_selector = logical(z_data_selector(z_data_selector(:,1) == current_track_ID,2:end));
    z_i = z(:,z_i_selector);
    R_i = R(:,:,z_i_selector);    
        
    if isempty(z_i) % || size(z_i,2) < 2
        continue;
    end % if    
        
    p_straight_ppd_i = nan(2,size(z_i,2));
    l_straight_ppd_i = nan(1,size(z_i,2));
    omega_r_straight_i = NaN(2,2,size(z_i,2));
    for j = 1:size(z_i,2)        
        % z_hat ___________________________________________________________
        l = dot(z_i(:,j)-p_0,t);
        p_straight_ppd_i(:,j) = p_0 + l*t;
        l_straight_ppd_i(:,j) = l;
        
        % omega ___________________________________________________________
        % omega_r_straight_i(:,:,j) = R_i(:,:,j)^-1;
        if det(R_i(:,:,j)) == 0
            omega_r_straight_i(:,:,j) = inf(2);
        else
            omega_r_straight_i(:,:,j) = invChol_mex(R_i(:,:,j));
        end % if      
    end % for j
    p_straight_ppd{i} = p_straight_ppd_i;
    l_straight_ppd{i} = l_straight_ppd_i;
    
%     %%% Test Plot
%     plot(z_i(1,:),z_i(2,:),'kx','MarkerSize',5); hold on;
%     plot(p_straight_ppd_i(1,:),p_straight_ppd_i(2,:),'r.','MarkerSize',5)
%     for test_i = 1:length(z_i)        
%         plot([z_i(1,test_i),p_straight_ppd_i(1,test_i)],[z_i(2,test_i),p_straight_ppd_i(2,test_i)],'b-');
%     end % for i   
%     axis equal
%     %%%^^^
    
    
    r_straight{i} = z_i - p_straight_ppd_i;
%     r_length{i} = (z_i(:,1)-z_i(:,end)) - (-L*t);
%     r_length{i} = [ norm(z_i(:,1)-p_0); norm(z_i(:,end)-(p_0+L*t)) ];
    r_length{i} = [ sqrt(norm(z_i(:,1)-p_0)^2 - norm(r_straight{i}(:,1))^2); sqrt(norm(z_i(:,end)-(p_0+L*t))^2 - norm(r_straight{i}(:,end))^2) ];
    omega_r_straight{i} = omega_r_straight_i;
    
%     if ~isempty(z_i)
%         r_length{i} = (z_i(:,1)-z_i(:,end)) - (-L*t);
%     else
%         r_length{i} = L*t;
%         %r_length{i} = [];
%     end % if

end % for i

end % function

