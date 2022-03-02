% 
% Plot data relating to localization
%   
%   Author: Hanno Winter
%   Date: 26-Mar-2021; Last revision: 17-May-2021

%% Position 

if(1)    
    
    plot_error_ellipses = 0; % enable or disable error-ellipse plots    
    plot_optim_data = 1; % plot optimized railway-map
    plot_ref_map_data = 1; % plot reference map
                    
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    if plot_ref_map_data
        prepareRefMapData
    end % if
    if plot_optim_data
        prepareOptimOutData    
    end % if

    % Plot ________________________________________________________________
    
    figure_name = 'Position';
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);
    clear h_plot    
    h_plot = gobjects(0);

    hold on; grid on;    
 
    h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg,sim_gnss_filter_input.Latitude_deg,'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
    h_plot(end+1) = plot(sim_imm.Longitude_deg,sim_imm.Latitude_deg,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_ekf.Longitude_deg,sim_ekf.Latitude_deg,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_tgc.Longitude_deg,sim_tgc.Latitude_deg,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    % h_plot(end+1) = plot(sim_tgc_map.Longitude_deg,sim_tgc_map.Latitude_deg,'r-','LineWidth',3,'MarkerSize',10,'DisplayName','TGC Tracks');
    h_plot(end+1) = plot(sim_tgc_map.Straights_Longitude_deg,sim_tgc_map.Straights_Latitude_deg,'-','Color',[1 0.5 0],'LineWidth',3,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
    h_plot(end+1) = plot(sim_tgc_map.Arcs_Longitude_deg,sim_tgc_map.Arcs_Latitude_deg,'-','Color',[1 1 1],'LineWidth',3,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');
    if plot_optim_data
        h_plot(end+1) = plot(optim_map_abs_position_longitude,optim_map_abs_position_latitude,'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');
        
%         test_idx1 = min(find(isnan(optim_map_abs_position_longitude))+1,length(optim_map_abs_position_longitude)-1);
%         test_idx2 = max(find(isnan(optim_map_abs_position_longitude))-1,1);
%         h_plot(end+1) = plot(optim_map_abs_position_longitude(test_idx1),optim_map_abs_position_latitude(test_idx1),'rx','LineWidth',2,'MarkerSize',10,'DisplayName','...');
%         h_plot(end+1) = plot(z_start_test_longitude,z_start_test_latitude,'r.','LineWidth',2,'MarkerSize',10,'DisplayName','...');
%         h_plot(end+1) = plot(optim_map_abs_position_longitude(test_idx2),optim_map_abs_position_latitude(test_idx2),'kx','LineWidth',2,'MarkerSize',10,'DisplayName','...');
%         h_plot(end+1) = plot(z_end_test_longitude,z_end_test_latitude,'k.','LineWidth',2,'MarkerSize',10,'DisplayName','...');
                
    end % if    
    if plot_ref_map_data && plot_ref_map
        h_plot(end+1) = plot(ref_track_map.Longitude_deg(ref_track_map_selector),ref_track_map.Latitude_deg(ref_track_map_selector),'y-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
    end % if      
    if plot_error_ellipses
        h_plot(end+1) = plot(sim_ekf_error_data.ee.Longitude_deg,sim_ekf_error_data.ee.Latitude_deg,'m-','LineWidth',1.5,'MarkerSize',10, ... 
            'DisplayName',sprintf('EKF %.2f Sigma-Error-Ellipse (appr. %.0f%%)',sim_ekf_error_data.ee.sigma_gain,sim_ekf_error_data.ee.confidence*100));
        h_plot(end+1) = plot(sim_imm_error_data.ee.Longitude_deg,sim_imm_error_data.ee.Latitude_deg,'g-','LineWidth',1.5,'MarkerSize',10, ... 
            'DisplayName',sprintf('IMM %.2f Sigma-Error-Ellipse (appr. %.0f%%)',sim_imm_error_data.ee.sigma_gain,sim_imm_error_data.ee.confidence*100));
        h_plot(end+1) = plot(sim_tgc_error_data.ee.Longitude_deg,sim_tgc_error_data.ee.Latitude_deg,'b-','LineWidth',1.5,'MarkerSize',10, ... 
            'DisplayName',sprintf('TGC %.2f Sigma-Error-Ellipse (appr. %.0f%%)',sim_tgc_error_data.ee.sigma_gain,sim_tgc_error_data.ee.confidence*100));
        h_plot(end+1) = plot(sim_gnss_filter_input_error_data.ee.Longitude_deg,sim_gnss_filter_input_error_data.ee.Latitude_deg,'k-','LineWidth',1.5,'MarkerSize',10, ... 
            'DisplayName',sprintf('GPS %.2f Sigma-Error-Ellipse (appr. %.0f%%)',sim_gnss_filter_input_error_data.ee.sigma_gain,sim_gnss_filter_input_error_data.ee.confidence*100));
    end % if
    
    if plot_ref_map_data && plot_ref_map && plot_optim_data
        xlim([min(optim_map_abs_position_longitude) max(optim_map_abs_position_longitude)])
        ylim([min(optim_map_abs_position_latitude) max(optim_map_abs_position_latitude)])
    end % if
    
    % Satellite plot
    % See: https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
    if exist('plot_google_map')
        plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
    end % if

    h_legend = legend(h_plot);
    set(h_legend,'Location','northeast')
    xlabel('longitude [deg]')
    ylabel('latitude [deg]')
    % axis equal  
    
    

end % if

%% Attitude

if(1)  
    
    % Plot-Calculations ___________________________________________________
    initLocalization
    prepareSimOutData
    prepareRefMapData    
            
    % Plot ________________________________________________________________
    figure_name = ['Attitude'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid on;

    clear h_plot    
    h_plot = gobjects(0);   
    ax1 = subplot(2,1,1); hold on; grid on;
    if plot_ref_map && ~strcmp(input_data_selector,'BS')
        h_plot(end+1) = plot(sim_imm.Time,imm_ref_track_map.Roll_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot(end+1) = plot(sim_attitude.Time,sim_attitude.Roll_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS');
    h_plot(end+1) = plot(sim_attitude_raw.Time,sim_attitude_raw.Roll_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS (no SD)');
    h_legend = legend(h_plot);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('roll [deg]')

    clear h_plot
    h_plot = gobjects(0); 
    ax2 = subplot(2,1,2); hold on; grid on;
    if plot_ref_map && ~strcmp(input_data_selector,'BS')
        h_plot(end+1) = plot(sim_imm.Time,imm_ref_track_map.Pitch_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot(end+1) = plot(sim_attitude.Time,sim_attitude.Pitch_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS');
    h_plot(end+1) = plot(sim_attitude_raw.Time,sim_attitude_raw.Pitch_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS (no SD)');
    h_legend = legend(h_plot);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('pitch [deg]')

    linkaxes([ax1,ax2],'x');
%     xlim([1 200])
        
end % if

%% INS Speeds

if(1)  
    
    % Plot-Calculations ___________________________________________________
	initLocalization
    prepareSimOutData
    
    gnss_v_ground = sqrt(sim_gnss.VelocityNorth_ms.^2+sim_gnss.VelocityEast_ms.^2);
    ins_v_ground = sqrt(sim_v_ned.VelocityNorth_ms.^2+sim_v_ned.VelocityEast_ms.^2);
        
    dir_corr_factor = abs(angdiff(deg2rad(sim_attitude.Heading_deg),deg2rad(sim_gnss.Heading_deg))) > pi/2;
    dir_corr_factor = dir_corr_factor * (-1);
    dir_corr_factor(dir_corr_factor==0) = 1;
    
    % Plot ________________________________________________________________
    figure_name = ['INS Speed Calculation'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid on;

    clear h_plot    
    h_plot = gobjects(0);   
    h_plot(end+1) = plot(sim_gnss.Time,dir_corr_factor.*gnss_v_ground,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
    h_plot(end+1) = plot(sim_ekf.Time,sim_ekf.VelocityVehicle_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(EKF)');
    h_plot(end+1) = plot(sim_imm.Time,sim_imm.VelocityVehicle_ms,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(IMM');
    h_plot(end+1) = plot(sim_a_eb_n.Time,sim_cs_v_ground_corrected,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(INS (with SD))');
    h_plot(end+1) = plot(sim_a_ib_b.Time,sim_cs_v_ground_uncorrected,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(INS (raw))');
    h_legend = legend(h_plot);
    set(h_legend,'Location','northwest')
    xlabel('time [s]')
    ylabel('v_{ground} [m/s]')
    
    % Plot ________________________________________________________________
    figure_name = ['INS Speed Calculation: Deltas'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid on;

    clear h_plot    
    h_plot = gobjects(0);   
    h_plot(end+1) = plot(sim_a_eb_n.Time,sim_cs_v_ground_corrected-dir_corr_factor.*gnss_v_ground,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','\Delta(INS with SD)');
    h_plot(end+1) = plot(sim_a_ib_b.Time,sim_cs_v_ground_uncorrected-dir_corr_factor.*gnss_v_ground,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','\Delta(INS raw)');
    h_legend = legend(h_plot);
    set(h_legend,'Location','northwest')
    xlabel('time [s]')
    ylabel('\Deltav_{ground} [m/s]')
        
end % if

%% Standstill detection

if(1)          
         
    % Plot-Calculations ___________________________________________________
    initLocalization
    prepareSimOutData
    
    gnss_speed = sqrt(sim_gnss.VelocityNorth_ms.^2+sim_gnss.VelocityEast_ms.^2+sim_gnss.VelocityDown_ms.^2);
    scaled_standstill_flag = sim_standstill_flag.standstill_flag*max(gnss_speed);
    
    % Plot ________________________________________________________________
    figure_name = ['Standstill Detection'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid on;

    clear h_plot    
    h_plot = gobjects(0);   
    ax1 = subplot(3,1,1); hold on; grid on;
    h_plot(end+1) = plot(sim_w_energy.Time,sim_w_energy.w_energy,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','E_w');
    h_plot(end+1) = plot(sim_w_energy.Time([1 end]),ones(1,2)*E_w_limit,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','E_w Limit');
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('E_{\omega} [rad^2/s^2]')

    clear h_plot
    h_plot = gobjects(0); 
    ax2 = subplot(3,1,2); hold on; grid on;
    h_plot(end+1) = plot(sim_a_e.Time,sim_a_e.a_e,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','a_e');
    h_plot(end+1) = plot(sim_a_e.Time([1 end]),ones(1,2)*a_e_limit,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','a_e Limit');
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('a_e [m/s^2]')

    clear h_plot
    h_plot = gobjects(0); 
    ax3 = subplot(3,1,3); hold on; grid on;
    h_plot(end+1) = plot(sim_gnss.Time,gnss_speed,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS Speed');
    h_plot(end+1) = plot(sim_ekf.Time,abs(sim_ekf.VelocityVehicle_ms),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(EKF)');
    h_plot(end+1) = plot(sim_imm.Time,abs(sim_imm.VelocityVehicle_ms),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(IMM)');
    h_plot(end+1) = plot(sim_standstill_flag.Time,scaled_standstill_flag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','(Scaled) Standstill Flag');
    h_legend = legend(h_plot);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('speed [m/s]')

    linkaxes([ax1,ax2,ax3],'x');
    % xlim(plot_time);
        
end % if

%% CDF: Uncertainty

if(1)  
    
    % Plot-Calculations ___________________________________________________
    initLocalization
    prepareSimOutData
    prepareRefMapData 
            
    % Plot ________________________________________________________________
    figure_name = ['CDF'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid on;
    
    title_str = sprintf([ ... 
                          'Localization Uncertainty\n', ... 
                          'conf.: ',num2str(error_ellipse_confidence),'%%\n', ... 
                          'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h\n', ... 
                          'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'\n', ... 
                          'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)) ... 
                       ]);
    sgtitle(title_str);
    
    clear h_plot
    h_plot = gobjects(0); 
    ax1 = subplot(1,3,1); hold on; grid on;
    %title(['Max (conf.: ',num2str(error_ellipse_confidence),'%, abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);    
    title('Max')
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.cdf_max.MaxError_m,sim_gnss_filter_input_error_data.cdf_max.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    h_plot(end+1) = plot(sim_ekf_error_data.cdf_max.MaxError_m,sim_ekf_error_data.cdf_max.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.cdf_max.MaxError_m,sim_imm_error_data.cdf_max.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.cdf_max.MaxError_m,sim_tgc_error_data.cdf_max.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel('error [m]')
    ylabel('availability')
%     xlim([0 10])
    ylim([0 1.2])
    
    clear h_plot
    h_plot = gobjects(0); 
    ax2 = subplot(1,3,2); hold on; grid on;
    %title(['AT (conf.: ',num2str(error_ellipse_confidence),'%, abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);
    title('AT')
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.cdf_at.MaxError_m,sim_gnss_filter_input_error_data.cdf_at.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    h_plot(end+1) = plot(sim_ekf_error_data.cdf_at.MaxError_m,sim_ekf_error_data.cdf_at.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.cdf_at.MaxError_m,sim_imm_error_data.cdf_at.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.cdf_at.MaxError_m,sim_tgc_error_data.cdf_at.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel('error [m]')
    ylabel('availability')
%     xlim([0 10])
    ylim([0 1.2])
    
    clear h_plot
    h_plot = gobjects(0); 
    ax3 = subplot(1,3,3); hold on; grid on; 
    %title(['CT (conf.: ',num2str(error_ellipse_confidence),', abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);
    title('CT')
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.cdf_ct.MaxError_m,sim_gnss_filter_input_error_data.cdf_ct.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    h_plot(end+1) = plot(sim_ekf_error_data.cdf_ct.MaxError_m,sim_ekf_error_data.cdf_ct.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.cdf_ct.MaxError_m,sim_imm_error_data.cdf_ct.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.cdf_ct.MaxError_m,sim_tgc_error_data.cdf_ct.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel('error [m]')
    ylabel('availability') 
%     xlim([0 10])
    ylim([0 1.2])
        
end % if

%% Perpendicular Deviation to Ref-Map

if(1)  
        
    % Plot-Calculations ___________________________________________________
	initLocalization
    prepareSimOutData
    prepareOptimOutData
    prepareRefMapData
        
    % Plot ________________________________________________________________
    if ~isempty(ref_track_map)    
        figure_name = ['Perpendicular Deviation to Ref-Map'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;
        % title(['Perpendicular deviation to reference-map (optim-map density: ',num2str(optim_map_density),'m, abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);
        title_str = sprintf([ ... 
                              'Perpendicular deviation to reference-map\n', ... 
                              'optim-map density: ',num2str(optim_map_density),'m\n', ... 
                              'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h\n', ... 
                              'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'\n', ... 
                              'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)) ... 
                           ]);
        title(title_str);
        
        clear h_plot    
        h_plot = gobjects(0);   
        h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_time,sim_gnss_filter_input_error_data.e_pp,'k.','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        h_plot(end+1) = plot(sim_ekf_error_data.e_time,sim_ekf_error_data.e_pp,'m.','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_plot(end+1) = plot(sim_imm_error_data.e_time,sim_imm_error_data.e_pp,'g.','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
        h_plot(end+1) = plot(sim_tgc_error_data.e_time,sim_tgc_error_data.e_pp,'b.','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
        %h_plot(end+1) = plot(optim_map_error_data.time,optim_map_error_data.e_pp,'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');

        h_legend = legend(h_plot);
        set(h_legend,'Location','southeast')
        xlabel('time [s]')
        ylabel('error [m]')
    end % if
end % if

%% CDF: Perpendicular Deviation to Ref-Map

if(1)  
        
    % Plot-Calculations ___________________________________________________
	initLocalization
    prepareSimOutData
    prepareOptimOutData
    prepareRefMapData
        
    % Plot ________________________________________________________________
    if ~isempty(ref_track_map)    
        figure_name = ['Perpendicular Deviation to Ref-Map CDF'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;
        
        title_str = sprintf([ ... 
                              'Perpendicular deviation to reference-map\n', ... 
                              'optim-map density: ',num2str(optim_map_density),'m\n', ... 
                              'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h\n', ... 
                              'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'\n', ... 
                              'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)) ... 
                           ]);
        title(title_str);

        clear h_plot    
        h_plot = gobjects(0);   
        h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_pp_cdf.MaxError_m,sim_gnss_filter_input_error_data.e_pp_cdf.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        h_plot(end+1) = plot(sim_ekf_error_data.e_pp_cdf.MaxError_m,sim_ekf_error_data.e_pp_cdf.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_plot(end+1) = plot(sim_imm_error_data.e_pp_cdf.MaxError_m,sim_imm_error_data.e_pp_cdf.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
        h_plot(end+1) = plot(sim_tgc_error_data.e_pp_cdf.MaxError_m,sim_tgc_error_data.e_pp_cdf.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
        h_plot(end+1) = plot(optim_map_error_data.e_pp_cdf.MaxError_m,optim_map_error_data.e_pp_cdf.Availability,'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');

        h_legend = legend(h_plot);
        set(h_legend,'Location','southeast')
        xlabel('error [m]')
        ylabel('availability')
        ylim([0 1.2])
                
    end % if
end % if

%%  Track-Length Deviation

if(1)  
        
    % Plot-Calculations ___________________________________________________
	initLocalization
    prepareSimOutData
    prepareOptimOutData
    prepareRefMapData
        
    % Plot ________________________________________________________________
    if ~isempty(ref_track_map)    
        figure_name = ['Track-Length Deviation'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;
        %title(['Track-length deviation to reference-map (optim-map density: ',num2str(optim_map_density),'m, abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);
        title_str = sprintf([ ... 
              'Track-length deviation to reference-map\n', ... 
              'optim-map density: ',num2str(optim_map_density),'m\n', ... 
              'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h\n', ... 
              'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'\n', ... 
              'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)) ... 
           ]);
        title(title_str);
        
        clear h_plot    
        h_plot = gobjects(0);   
        h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_time,sim_gnss_filter_input_error_data.e_tl,'k.','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        h_plot(end+1) = plot(sim_ekf_error_data.e_time,sim_ekf_error_data.e_tl,'m.','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_plot(end+1) = plot(sim_imm_error_data.e_time,sim_imm_error_data.e_tl,'g.','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
        h_plot(end+1) = plot(sim_tgc_error_data.e_time,sim_tgc_error_data.e_tl,'b.','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
%         h_plot(end+1) = plot(optim_map_error_data.e_tl_cdf.MaxError_m,optim_map_error_data.e_tl_cdf.Availability,'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');

        h_legend = legend(h_plot);
        set(h_legend,'Location','southeast')
        xlabel('time [s]')
        ylabel('error [m]')
        
    end % if
end % if

%% CDF: Track-Length Deviation

if(1)  
        
    % Plot-Calculations ___________________________________________________
	initLocalization
    prepareSimOutData
    prepareOptimOutData
    prepareRefMapData
        
    % Plot ________________________________________________________________
    if ~isempty(ref_track_map)    
        figure_name = ['Track-Length Deviation CDF'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;
        % title(['Track-length deviation to reference-map (optim-map density: ',num2str(optim_map_density),'m, abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);

        title_str = sprintf([ ... 
                      'Track-length deviation to reference-map\n', ... 
                      'optim-map density: ',num2str(optim_map_density),'m\n', ... 
                      'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h\n', ... 
                      'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'\n', ... 
                      'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)) ... 
                   ]);
        title(title_str);
        
        clear h_plot    
        h_plot = gobjects(0);   
        h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_tl_cdf.MaxError_m,sim_gnss_filter_input_error_data.e_tl_cdf.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        h_plot(end+1) = plot(sim_ekf_error_data.e_tl_cdf.MaxError_m,sim_ekf_error_data.e_tl_cdf.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
        h_plot(end+1) = plot(sim_imm_error_data.e_tl_cdf.MaxError_m,sim_imm_error_data.e_tl_cdf.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
        h_plot(end+1) = plot(sim_tgc_error_data.e_tl_cdf.MaxError_m,sim_tgc_error_data.e_tl_cdf.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
        h_plot(end+1) = plot(optim_map_error_data.e_tl_cdf.MaxError_m,optim_map_error_data.e_tl_cdf.Availability,'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');

        h_legend = legend(h_plot);
        set(h_legend,'Location','southeast')
        xlabel('error [m]')
        ylabel('availability')
        ylim([0 1.2])
        
        % Output ______________________________________________________________

        fprintf('\n### Max Errors ###\n\n') 
        fprintf('GPS: %.2f\n',max(sim_gnss_filter_input_error_data.e_tl));
        fprintf('EKF: %.2f\n',max(sim_ekf_error_data.e_tl));
        fprintf('IMM: %.2f\n',max(sim_imm_error_data.e_tl));
        fprintf('TGC: %.2f\n',max(sim_tgc_error_data.e_tl));
        fprintf('Map: %.2f\n',max(optim_map_error_data.e_tl));
        
    end % if
end % if

%% Helper functions

function [h_tracks, h_map] = plotXyTrackMap(xy_track_map)
% [h_tracks, h_map] = plotXyTrackMap(xy_tack_map)
%

% Init ____________________________________________________________________

tracks = unique(xy_track_map(:,1).track_num);
h_tracks = gobjects(length(tracks),1);
map_data = table();

% Calculations ____________________________________________________________

for i = 1:length(tracks)
    
    track_i = tracks(i);    
    track_i_selector = (xy_track_map(:,1).track_num == track_i);    
    track_i_data = xy_track_map(track_i_selector,:);
    
    track_i_map_data = track_i_data;
    track_i_map_data{end+1,:} = nan;
    map_data = [map_data; track_i_map_data];
    
%     if verLessThan('matlab', '9.5')     
        h_tracks(i) = plot(track_i_data.lon,track_i_data.lat,'DisplayName',sprintf('Track: %i',track_i)); hold on;
%     else
%         h_tracks(i) = geoplot(track_i_data.lat,track_i_data.lon,'DisplayName',sprintf('Track: %i',track_i)); hold on;
%         geobasemap('topographic')        
%     end % if
    
end % for i

% Plot whole map at once
% if verLessThan('matlab', '9.5') 
    h_map = plot(map_data.lon,map_data.lat);
% else
%     h_map = geoplot(map_data.lat,map_data.lon);
% end % if

end % function