% 
% Plot data to be shown in thesis
%   
%   Author: Hanno Winter
%   Date: 03-Jul-2021; Last revision: 15-Dec-2021

%% INS Speed Error

if(1)
    
    initLocalization
    prepareSimOutData
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            plot_start_time = 0; 
            plot_end_time = 180;
        case {'BS'}
            plot_start_time = 0; 
            plot_end_time = 180;
        case {'NT'}
            plot_start_time = 0; 
            plot_end_time = 180;
    end % switch
    
    ds_rate = 1; % in sec
                              
    % Plot-Calculations ___________________________________________________  
    
    ds_sim_data = max(ceil(ds_rate/imu_sample_time),1);
    ds_filter_data = max(ceil(ds_rate/imm_sample_time),1);
    ds_gnss_data = max(ceil(ds_rate/gnss_sample_time),1);
    
    if use_time_frame
        plot_selector = (sim_a_ib_b.Time >= plot_start_time) & (sim_a_ib_b.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_a_ib_b.Time));
    end % if
    plot_indices = find(plot_selector);
    
    x_limits = [min(sim_a_ib_b.Time(plot_selector)) max(sim_a_ib_b.Time(plot_selector))];
    
    standstill_str = sprintf('%d',sim_standstill_flag.standstill_flag);
    standstill_start_idx = strfind(standstill_str,'01')'+1;
    if sim_standstill_flag.standstill_flag(1) == true
        standstill_start_idx = [1;standstill_start_idx];
    end % if
    standstill_end_idx = strfind(standstill_str,'10')';
    if sim_standstill_flag.standstill_flag(end) == true
        standstill_end_idx = [standstill_end_idx;length(sim_standstill_flag.standstill_flag)];
    end % if
    standstill_start_times = sim_a_ib_b.Time(standstill_start_idx)-x_limits(1);
    standstill_end_times = sim_a_ib_b.Time(standstill_end_idx)-x_limits(1);
    standstill_border_times = [standstill_start_times,standstill_end_times];
                
    dir_corr_factor = abs(angdiff(deg2rad(sim_attitude.Heading_deg),deg2rad(sim_gnss.Heading_deg))) > pi/2;
    dir_corr_factor = dir_corr_factor * (-1);
    dir_corr_factor(dir_corr_factor==0) = 1;
    
    gnss_raw_data_rs = resample(gnss_sim_data,gnss_sim_data.Time(1):imu_sample_time:gnss_sim_data.Time(end),'zoh');
    gnss_v_ground = sqrt(gnss_raw_data_rs.Data(:,4).^2+gnss_raw_data_rs.Data(:,5).^2);  
%     gnss_v_ground = sqrt(sim_gnss.VelocityNorth_ms.^2+sim_gnss.VelocityEast_ms.^2);
    ins_v_ground = sqrt(sim_v_ned.VelocityNorth_ms.^2+sim_v_ned.VelocityEast_ms.^2);
    
%     e_v_ins = abs(abs(sim_cs_v_ground_corrected)-gnss_v_ground) ./ max(gnss_v_ground);
    e_v_ins = abs(sim_cs_v_ground_xy_corrected-gnss_v_ground(1:length(sim_cs_v_ground_xy_corrected)));
%     e_v_raw = abs(abs(sim_cs_v_ground_uncorrected)-gnss_v_ground) ./ max(gnss_v_ground);
    e_v_raw = abs(sim_cs_v_ground_uncorrected-dir_corr_factor.*gnss_v_ground(1:length(sim_cs_v_ground_xy_corrected)));    
    
%     time_ds = (0:gnss_sample_time:sim_a_ib_b.Time(end))';
%     e_v_ins_ds = interp1(sim_a_eb_n.Time,e_v_ins,time_ds);
%     e_v_raw_ds = interp1(sim_a_ib_b.Time,e_v_raw,time_ds);
            
    % Plot ________________________________________________________________
    figure_name = ['INS Speed Error'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid off;

    clear h_plot    
    h_plot = gobjects(0);
    %h_plot(end+1) = plot(sim_gnss.Time-x_limits(1),dir_corr_factor.*gnss_v_ground,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
    h_plot(end+1) = plot(sim_a_eb_n.Time(1:ds_sim_data:end)-x_limits(1),e_v_ins(1:ds_sim_data:end)*3.6,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','mit Lagekorrektur');
    h_plot(end+1) = plot(sim_a_ib_b.Time(1:ds_sim_data:end)-x_limits(1),e_v_raw(1:ds_sim_data:end)*3.6,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','ohne Lagekorrektur');
    axis tight
    
    if exist('x_limits','var')
        xlim([0 x_limits(2)-x_limits(1)])
%         y_limits = [min([e_v_ins(plot_selector);e_v_raw(plot_selector)]) max([e_v_ins(plot_selector);e_v_raw(plot_selector)])];
%         ylim(y_limits);
    end % if
    
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','northwest')
    xlabel('time [s]')
    ylabel('v_{GPS} [km/h]')
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/v_estim_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/v_estim_',input_data_selector], ...
                    'relativeDataPath',['figs/data/v_estim_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
        
end % if

%% Attitude

if(1)
	
    initLocalization
    prepareSimOutData
    prepareRefMapData
    
    use_time_frame = 0;
    switch input_data_selector
        case {'C'}
            plot_start_time = 0; 
            plot_end_time = 885;
        case {'BS'}
            plot_start_time = 0; % 734; 
            plot_end_time = 330; % 1040;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    ds_rate = 0.1; % in sec
                           
    % Plot-Calculations ___________________________________________________  
    
    ds_sim_data = max(ceil(ds_rate/imu_sample_time),1);
    ds_filter_data = max(ceil(ds_rate/imm_sample_time),1);
    ds_gnss_data = max(ceil(ds_rate/gnss_sample_time),1);
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
    end % if
    plot_indices = find(plot_selector);
    
    x_limits = [min(sim_imm.Time(plot_selector)) max(sim_imm.Time(plot_selector))];
    
    standstill_str = sprintf('%d',sim_standstill_flag.standstill_flag);
    standstill_start_idx = strfind(standstill_str,'01')'+1;
    if sim_standstill_flag.standstill_flag(1) == true
        standstill_start_idx = [1;standstill_start_idx];
    end % if
    standstill_end_idx = strfind(standstill_str,'10')';
    if sim_standstill_flag.standstill_flag(end) == true
        standstill_end_idx = [standstill_end_idx;length(sim_standstill_flag.standstill_flag)];
    end % if
    standstill_start_times = sim_a_ib_b.Time(standstill_start_idx)-x_limits(1);
    standstill_end_times = sim_a_ib_b.Time(standstill_end_idx)-x_limits(1);
    standstill_border_times = [standstill_start_times,standstill_end_times];
    
    gnss_v_ground = sqrt(sim_gnss.VelocityNorth_ms.^2+sim_gnss.VelocityEast_ms.^2);
    
    time_ds = (0:gnss_sample_time:sim_gnss.Time(end))';
    gnss_v_ground_ds = interp1(sim_gnss.Time,gnss_v_ground,time_ds);
    
    % Plot ________________________________________________________________
    figure_name = ['Attitude'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid off;
   
    h_plot_ax1 = gobjects(0);   
    ax1 = subplot(3,1,1); hold on; grid off;
    if plot_ref_map && ~strcmp(input_data_selector,'BS')
        h_plot_ax1(end+1) = plot(sim_imm.Time(1:ds_filter_data:end)-x_limits(1),imm_ref_track_map.Roll_deg(1:ds_filter_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot_ax1(end+1) = plot(sim_attitude.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude.Roll_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS');
    h_plot_ax1(end+1) = plot(sim_attitude_raw.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude_raw.Roll_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS (no SD)');
    axis tight    
   
    h_plot_ax2 = gobjects(0); 
    ax2 = subplot(3,1,2); hold on; grid off;
    if plot_ref_map && ~strcmp(input_data_selector,'BS')
        h_plot_ax2(end+1) = plot(sim_imm.Time(1:ds_filter_data:end)-x_limits(1),imm_ref_track_map.Pitch_deg(1:ds_filter_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot_ax2(end+1) = plot(sim_attitude.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude.Pitch_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS');
    h_plot_ax2(end+1) = plot(sim_attitude_raw.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude_raw.Pitch_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS (no SD)');
    axis tight
        
    h_plot_ax3 = gobjects(0); 
    ax3 = subplot(3,1,3); hold on; grid off;
    %h_plot_ax3(end+1) = plot(sim_imm.Time-x_limits(1),sim_imm.VelocityVehicle_ms*3.6,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot_ax3(end+1) = plot(sim_gnss.Time(1:ds_sim_data:end)-x_limits(1),gnss_v_ground(1:ds_sim_data:end)*3.6,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    axis tight
   
    linkaxes([ax1,ax2,ax3],'x');
    
    if exist('x_limits','var')
        %axis equal
        xlim([0 x_limits(2)-x_limits(1)])
        %ylim(y_limits)
    end % if
    
    axes(ax1)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    h_legend = legend(h_plot_ax1);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('roll [deg]')
    
    axes(ax2)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    h_legend = legend(h_plot_ax2);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('pitch [deg]')
    
    axes(ax3)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');    
    h_legend = legend(h_plot_ax3);
    set(h_legend,'Location','southeast')
    xlabel('Zeit in s')
    ylabel('v in km/h')
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/attitude_estim_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/attitude_estim_',input_data_selector], ...
                    'relativeDataPath',['figs/data/attitude_estim_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
        
end % if

%% Initial Attitude Offset

if(1)
	
    initLocalization
    prepareSimOutData
    prepareRefMapData
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            plot_start_time = 0; 
            plot_end_time = 3;
        case {'BS'}
            plot_start_time = 10; 
            plot_end_time = 100;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    ds_rate = 0.01; % in sec
                           
    % Plot-Calculations ___________________________________________________  
    
    ds_sim_data = max(ceil(ds_rate/imu_sample_time),1);
    ds_filter_data = max(ceil(ds_rate/imm_sample_time),1);
    ds_gnss_data = max(ceil(ds_rate/gnss_sample_time),1);
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
    end % if
    plot_indices = find(plot_selector);
    
    x_limits = [min(sim_imm.Time(plot_selector)) max(sim_imm.Time(plot_selector))];
    
    standstill_str = sprintf('%d',sim_standstill_flag.standstill_flag);
    standstill_start_idx = strfind(standstill_str,'01')'+1;
    if sim_standstill_flag.standstill_flag(1) == true
        standstill_start_idx = [1;standstill_start_idx];
    end % if
    standstill_end_idx = strfind(standstill_str,'10')';
    if sim_standstill_flag.standstill_flag(end) == true
        standstill_end_idx = [standstill_end_idx;length(sim_standstill_flag.standstill_flag)];
    end % if
    standstill_start_times = sim_a_ib_b.Time(standstill_start_idx)-x_limits(1);
    standstill_end_times = sim_a_ib_b.Time(standstill_end_idx)-x_limits(1);
    standstill_border_times = [standstill_start_times,standstill_end_times];
            
    % Plot ________________________________________________________________
    figure_name = ['Attitude: with initial offset'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid off;
   
    h_plot_ax1 = gobjects(0);   
    ax1 = subplot(2,1,1); hold on; grid off;
    if plot_ref_map && ~strcmp(input_data_selector,'BS')
        h_plot_ax1(end+1) = plot(sim_imm.Time-x_limits(1),imm_ref_track_map.Roll_deg,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot_ax1(end+1) = plot(sim_attitude.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude.Roll_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS');
    h_plot_ax1(end+1) = plot(sim_attitude_raw.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude_raw.Roll_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS (no SD)');
    axis tight    
   
    h_plot_ax2 = gobjects(0); 
    ax2 = subplot(2,1,2); hold on; grid off;
    if plot_ref_map && ~strcmp(input_data_selector,'BS')
        h_plot_ax2(end+1) = plot(sim_imm.Time(1:ds_filter_data:end)-x_limits(1),imm_ref_track_map.Pitch_deg(1:ds_filter_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot_ax2(end+1) = plot(sim_attitude.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude.Pitch_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS');
    h_plot_ax2(end+1) = plot(sim_attitude_raw.Time(1:ds_sim_data:end)-x_limits(1),sim_attitude_raw.Pitch_deg(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','INS (no SD)');
    axis tight
   
    linkaxes([ax1,ax2],'x');
    
    if exist('x_limits','var')
        %axis equal
        xlim([0 x_limits(2)-x_limits(1)])
        %ylim(y_limits)
    end % if
    
    axes(ax1)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
%     vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    h_legend = legend(h_plot_ax1);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('roll [deg]')
    
    axes(ax2)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
%     vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    h_legend = legend(h_plot_ax2);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('pitch [deg]')
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/attitude_offset_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/attitude_offset_',input_data_selector], ...
                    'relativeDataPath',['figs/data/attitude_offset_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
        
end % if

%% Standstill detection

if(1)          
    
    initLocalization
    prepareSimOutData
    
    use_time_frame = 0;
    switch input_data_selector
        case {'C'}
            plot_start_time = 10; 
            plot_end_time = 120;
        case {'BS'}
            plot_start_time = 3; 
            plot_end_time = 300;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    ds_rate = 0.1; % in sec
            
    % Plot-Calculations ___________________________________________________
    
    ds_sim_data = max(ceil(ds_rate/imu_sample_time),1);
    ds_filter_data = max(ceil(ds_rate/imm_sample_time),1);
    ds_gnss_data = max(ceil(ds_rate/gnss_sample_time),1);
    
    if use_time_frame
        plot_selector = (sim_a_ib_b.Time >= plot_start_time) & (sim_a_ib_b.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_a_ib_b.Time));
    end % if
    plot_indices = find(plot_selector);
    
    x_limits = [min(sim_a_ib_b.Time(plot_selector)) max(sim_a_ib_b.Time(plot_selector))];
    
    standstill_str = sprintf('%d',sim_standstill_flag.standstill_flag);
    standstill_start_idx = strfind(standstill_str,'01')'+1;
    if sim_standstill_flag.standstill_flag(1) == true
        standstill_start_idx = [1;standstill_start_idx];
    end % if
    standstill_end_idx = strfind(standstill_str,'10')';
    if sim_standstill_flag.standstill_flag(end) == true
        standstill_end_idx = [standstill_end_idx;length(sim_standstill_flag.standstill_flag)];
    end % if
    standstill_start_times = sim_a_ib_b.Time(standstill_start_idx)-x_limits(1);
    standstill_end_times = sim_a_ib_b.Time(standstill_end_idx)-x_limits(1);
    standstill_border_times = [standstill_start_times,standstill_end_times];
    
    %gnss_speed = sqrt(sim_gnss.VelocityNorth_ms.^2+sim_gnss.VelocityEast_ms.^2);
    gnss_speed = sim_gnss_filter_input.VelocityGround_ms;
%     scaled_standstill_flag = sim_standstill_flag.standstill_flag*max(gnss_speed);
    
    % Plot ________________________________________________________________
    figure_name = ['Standstill Detection'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid on;
 
    h_plot_ax1 = gobjects(0);   
    ax1 = subplot(3,1,1); hold on; grid on;
    h_plot_ax1(end+1) = plot(sim_w_energy.Time(1:ds_sim_data:end)-x_limits(1),sim_w_energy.w_energy(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','E_w');
    h_plot_ax1(end+1) = plot(sim_w_energy.Time([1 end])-x_limits(1),ones(1,2)*E_w_limit,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','E_w_{TH}');
    axis tight    

    h_plot_ax2 = gobjects(0); 
    ax2 = subplot(3,1,2); hold on; grid on;
    h_plot_ax2(end+1) = plot(sim_a_e.Time(1:ds_sim_data:end)-x_limits(1),sim_a_e.a_e(1:ds_sim_data:end),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','a_e');
    h_plot_ax2(end+1) = plot(sim_a_e.Time([1 end])-x_limits(1),ones(1,2)*a_e_limit,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','a_e_{TH}');
    axis tight
    
    h_plot_ax3 = gobjects(0); 
    ax3 = subplot(3,1,3); hold on; grid on;
    h_plot_ax3(end+1) = plot(sim_gnss_filter_input.Time(1:ds_gnss_data:end)-x_limits(1),gnss_speed(1:ds_gnss_data:end)*3.6,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    %h_plot_ax3(end+1) = plot(sim_ekf.Time-x_limits(1),abs(sim_ekf.VelocityVehicle_ms),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(EKF)');
    %h_plot_ax3(end+1) = plot(sim_imm.Time-x_limits(1),abs(sim_imm.VelocityVehicle_ms),'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','abs(IMM)');
    %h_plot_ax3(end+1) = plot(sim_standstill_flag.Time-x_limits(1),scaled_standstill_flag,'-','LineWidth',1.5,'MarkerSize',10,'DisplayName','(Scaled) Standstill Flag');
    axis tight
    
    linkaxes([ax1,ax2,ax3],'x');
    % xlim(plot_time);
    
    if exist('x_limits','var')
        xlim([0 x_limits(2)-x_limits(1)])
    end % if    
    
    axes(ax1)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [0-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    h_legend = legend(h_plot_ax1);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('E_{\omega} [rad^2/s^2]')
    
    axes(ax2)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [0-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');
    h_legend = legend(h_plot_ax2);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('a_e [m/s^2]')
    
    axes(ax3)
    %y_limits = ylim; delta_y_lim = diff(y_limits);
    %y_limits_new = [0-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim([0 20])
    vfill(standstill_border_times,'b','facealpha',.2,'edgecolor','none');    
    h_legend = legend(h_plot_ax3);
    set(h_legend,'Location','southwest')
    xlabel('time [s]')
    ylabel('speed [km/h]') 
    
    clear h_plot_ax1 h_plot_ax2 h_plot_ax3
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/standstill_det_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/standstill_det_',input_data_selector], ...
                    'relativeDataPath',['figs/data/standstill_det_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
        
end % if

%% Position 

if(1)
    
    plot_error_ellipses = 0; % enable or disable error-ellipse plots    
    plot_optim_data = 1; % plot optimized railway-map
    plot_ref_map_data = 1; % plot reference map
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            plot_start_time = 1380; 
            plot_end_time = 1900;
        case {'BS'}
            plot_start_time = 10; 
            plot_end_time = 100;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
    end % if
    plot_indices = find(plot_selector);    
            
    x_limits = [min(sim_imm.Longitude_deg(plot_selector)) max(sim_imm.Longitude_deg(plot_selector))];
    y_limits = [min(sim_imm.Latitude_deg(plot_selector)) max(sim_imm.Latitude_deg(plot_selector))];
                    
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
    
    title_str = sprintf(['Positions\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
                       ]);
    sgtitle(title_str);
    
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
    
    if exist('x_limits','var')
        %axis equal
        xlim(x_limits)
        ylim(y_limits)
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

%% Geometry Identification

if(1)
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
%             plot_start_time = 5; 
%             plot_end_time = 800;
            plot_start_time = 1000; 
            plot_end_time = 1380;
        case {'BS'}
            plot_start_time = 750; 
            plot_end_time = 900;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_indices = find(plot_selector);
    
    sim_map_track_element_selector = false(size(sim_map.track_maps.ID));
    for i = 1:length(sim_map.track_maps.ID)
        track_id_i = sim_map.track_maps.ID(i);
        optimization_data_selector_idx = find(track_id_i==optimization_imm_data_selector(:,1));
                
        if any(optimization_imm_data_selector(optimization_data_selector_idx,2:end) & plot_selector(:)')
            sim_map_track_element_selector(i) = true;
        end % if
    end % for i
%     
%     sim_map.track_maps(sim_map_track_element_selector,:)
%     [sim_map.track_start_points(sim_map_track_element_selector,:).ID sim_map.track_start_points(sim_map_track_element_selector,:).x_0]
%     [sim_map.track_start_points(sim_map_track_element_selector,:).ID sim_map.track_start_points(sim_map_track_element_selector,:).y_0]
%     [sim_map.track_start_points(sim_map_track_element_selector,:).ID 90-sim_map.track_start_points(sim_map_track_element_selector,:).phi_0]
    
    tgc_map_track_element_selector = false(size(sim_tgc_map.track_id));
    for i = 1:length(sim_tgc_map.track_id)
        track_id_i = sim_tgc_map.track_id(i);
        optimization_data_selector_idx = find(track_id_i==optimization_imm_data_selector(:,1));
                
        if any(optimization_imm_data_selector(optimization_data_selector_idx,2:end) & plot_selector(:)')
            tgc_map_track_element_selector(i) = true;
        end % if
    end % for i
    
    straight_selector_temp = ismember(sim_map_straights.track_start_points.ID,sim_tgc_map.track_id(tgc_map_track_element_selector)); 
    sim_map_straights_temp.topology = sim_map_straights.topology(straight_selector_temp,straight_selector_temp);
    sim_map_straights_temp.track_start_points = sim_map_straights.track_start_points(straight_selector_temp,:);
    sim_map_straights_temp.track_maps = sim_map_straights.track_maps(straight_selector_temp,:);    
    [~,~,sim_map_straights_temp_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_straights_temp,sim_map_density);
    sim_map_straights_temp_abs_position_utm = p_0_utm(1:2) + sim_map_straights_temp_abs_position_utm;
    [sim_map_straights_temp_abs_position_lat,sim_map_straights_temp_abs_position_lon] = ... 
        utm2ll(sim_map_straights_temp_abs_position_utm(1,:),sim_map_straights_temp_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
    arc_selector_temp = ismember(sim_map_arcs.track_start_points.ID,sim_tgc_map.track_id(tgc_map_track_element_selector)); 
    sim_map_arc_temp.topology = sim_map_arcs.topology(arc_selector_temp,arc_selector_temp);
    sim_map_arc_temp.track_start_points = sim_map_arcs.track_start_points(arc_selector_temp,:);
    sim_map_arc_temp.track_maps = sim_map_arcs.track_maps(arc_selector_temp,:);    
    [~,~,sim_map_arc_temp_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_arc_temp,sim_map_density);
    sim_map_arc_temp_abs_position_utm = p_0_utm(1:2) + sim_map_arc_temp_abs_position_utm;
    [sim_map_arcs_temp_abs_position_lat,sim_map_arcs_temp_abs_position_lon] = ... 
        utm2ll(sim_map_arc_temp_abs_position_utm(1,:),sim_map_arc_temp_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
    [~,mu_max_idx] = max(sim_imm.ModelProbabilities,[],2);
    mu_geometry_selector = [(mu_max_idx == 1), (mu_max_idx == 2), (mu_max_idx == 3)];
    mu_geometry_selector = mu_geometry_selector(plot_selector,:);
    mu_straight_selector = (mu_max_idx == 1);
    mu_arc_selector = (mu_max_idx == 2);
    mu_unkown_selector = (mu_max_idx == 3);
    
    %mu_straight_time = sim_imm.Time(mu_straight_selector & plot_selector);
    %mu_arc_time = sim_imm.Time(mu_arc_selector & plot_selector);
    mu_straight_selector = mu_straight_selector(plot_selector);    
    mu_arc_selector = mu_arc_selector(plot_selector);
    mu_unkown_selector = mu_unkown_selector(plot_selector);
    
    mu_geometry_border_times = {};
    for i = 1:size(mu_geometry_selector,2)
        mu_geometry_str = sprintf('%d',mu_geometry_selector(:,i));
        mu_geometry_start_idx = strfind(mu_geometry_str,'01')'+1;
        if mu_geometry_selector(1,i) == true
            mu_geometry_start_idx = [1;mu_geometry_start_idx];
        end % if
        mu_geometry_end_idx = strfind(mu_geometry_str,'10')';
        if mu_geometry_selector(end,i) == true
            mu_geometry_end_idx = [mu_geometry_end_idx;length(mu_geometry_selector(:,i))];
        end % if
        mu_geometry_start_times = sim_imm.Time(mu_geometry_start_idx)-sim_imm.Time(1);
        mu_geometry_end_times = sim_imm.Time(mu_geometry_end_idx)-sim_imm.Time(1);
        mu_geometry_border_times{i} = [mu_geometry_start_times,mu_geometry_end_times];
    end % for i                
    
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Plot ________________________________________________________________
    
    figure_name = ['Geometry Identification I (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);      
    clear h_plot    
        
    clear h_plot    
    h_plot = gobjects(0); hold on; grid on;      
    if plot_ref_map
        h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector),sim_gnss_filter_input.Latitude_deg(plot_gnss_selector),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
        axis tight
        x_limits = h_plot(end).Parent.XLim;
        y_limits = h_plot(end).Parent.YLim;
        %h_plot(end+1) = plot(imm_ref_track_map.Longitude_deg(plot_selector),imm_ref_track_map.Latitude_deg(plot_selector),'k-','LineWidth',5,'MarkerSize',10,'DisplayName','Track-Map');  
        h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'k-','LineWidth',5,'MarkerSize',10,'DisplayName','Track-Map');  
        % h_plot(end+1) = plot(imm_ref_track_map.Longitude_deg(plot_indices(1)),imm_ref_track_map.Latitude_deg(plot_indices(1)),'kx','MarkerSize',10,'DisplayName','Start');  
        % h_plot(end+1) = plot(imm_ref_track_map.Longitude_deg(plot_indices(end)),imm_ref_track_map.Latitude_deg(plot_indices(end)),'ko','MarkerSize',10,'DisplayName','Stop');  
    end % if
    h_plot(end+1) = plot(sim_map_straights_temp_abs_position_lon,sim_map_straights_temp_abs_position_lat,'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
    h_plot(end+1) = plot(sim_map_arcs_temp_abs_position_lon,sim_map_arcs_temp_abs_position_lat,'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');
    if exist('x_limits','var')
        %axis equal
        xlim(x_limits)
        ylim(y_limits)
    end % if    
    
    % Annotation
    nan_idx = find(isnan(sim_map_straights_temp_abs_position_lon));
    mean_start_interval_idx = [1,nan_idx(1:end-1)+1];
    mean_end_interval_idx = nan_idx-1;
    mean_intervals = [mean_start_interval_idx',mean_end_interval_idx'];
    mean_straight_pos_lon = nan(length(mean_intervals),1);
    mean_straight_pos_lat = nan(length(mean_intervals),1);
    for i = 1:length(mean_intervals)
        mean_interval_idx = [mean_intervals(i,1):mean_intervals(i,2)];
        mean_straight_pos_lon(i) = mean(sim_map_straights_temp_abs_position_lon(mean_interval_idx));
        mean_straight_pos_lat(i) = mean(sim_map_straights_temp_abs_position_lat(mean_interval_idx));
    end % for i   
    text(mean_straight_pos_lon,mean_straight_pos_lat,num2str(sim_map_straights_temp.track_maps.ID))
    
    nan_idx = find(isnan(sim_map_arcs_temp_abs_position_lon));
    mean_start_interval_idx = [1,nan_idx(1:end-1)+1];
    mean_end_interval_idx = nan_idx-1;
    mean_intervals = [mean_start_interval_idx',mean_end_interval_idx'];
    mean_unkown_pos_lon = nan(length(mean_intervals),1);
    mean_unkown_pos_lat = nan(length(mean_intervals),1);
    for i = 1:length(mean_intervals)
        mean_interval_idx = [mean_intervals(i,1):mean_intervals(i,2)];
        mean_unkown_pos_lon(i) = mean(sim_map_arcs_temp_abs_position_lon(mean_interval_idx));
        mean_unkown_pos_lat(i) = mean(sim_map_arcs_temp_abs_position_lat(mean_interval_idx));
    end % for i   
    text(mean_unkown_pos_lon,mean_unkown_pos_lat,num2str(sim_map_arc_temp.track_maps.ID))
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','northeast')
    xlabel(['Longitude in deg'])
    ylabel(['Latitude in deg'])
%     xlabel(['L [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
%     ylabel(['y [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
    
    % Satellite plot
    % See: https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
    if exist('plot_google_map')
        plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
    end % if
    
    sgtitle(sprintf([ ... 
            'Geometry Identification I\n', ... 
            '(Data Set: ',input_data_selector,'; ', ... 
            '   Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
            '     End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
          ]));
      
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/geom_det_',input_data_selector,'_positions.tex'], ... 
                    'dataPath',['plots/paper_plots/data/geom_det_',input_data_selector,'_positions'], ...
                    'relativeDataPath',['figs/data/geom_det_',input_data_selector,'_positions'], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
    % _____________________________________________________________________
    
    figure_name = ['Geometry Identification II (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);    
    clear h_plot 
    
    clear h_plot    
    h_plot = gobjects(0);
    h_background = gobjects(0);
    ax1 = subplot(3,1,1); hold on; grid on;
    h_plot(end+1) = plot(sim_imu_filter_input.Time(plot_selector)-sim_imm.Time(plot_indices(1)),sim_imu_filter_input.TurnRateZ_degs(plot_selector),'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','w_z');
    axis tight
    ylim([-0.07 0.07])
    vfill(mu_geometry_border_times{1},'b','facealpha',.2,'edgecolor','none');
    vfill(mu_geometry_border_times{2},'g','facealpha',.2,'edgecolor','none');
    vfill(mu_geometry_border_times{3},'r','facealpha',.2,'edgecolor','none');   
    
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','northeast')
    xlabel('t [s]')
    ylabel('w_z [deg/s]')
    %axis tight
    ylim([-0.07 0.07])
    
    clear h_plot    
    h_plot = gobjects(0);   
    ax2 = subplot(3,1,2); hold on; grid on; 
    h_plot(end+1) = plot(sim_imm.Time(plot_selector)-sim_imm.Time(plot_indices(1)),sim_imm.ModelProbabilities(plot_selector,1),'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','straight');
    h_plot(end+1) = plot(sim_imm.Time(plot_selector)-sim_imm.Time(plot_indices(1)),sim_imm.ModelProbabilities(plot_selector,2),'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','circular arc');
    h_plot(end+1) = plot(sim_imm.Time(plot_selector)-sim_imm.Time(plot_indices(1)),sim_imm.ModelProbabilities(plot_selector,3),'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','unkown');
    h_legend = legend(h_plot);
    set(h_legend,'Location','northeast')
    xlabel('t [s]')
    ylabel('probability [-]')
    axis tight
    
    clear h_plot    
    h_plot = gobjects(0);   
    ax3 = subplot(3,1,3); hold on; grid on; 
    h_plot(end+1) = plot(sim_imu_filter_input.Time(plot_selector)-sim_imm.Time(plot_indices(1)),abs(sim_imm.VelocityVehicle_ms(plot_selector))*3.6,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','v_{vehicle}');
    v_vehicle_kmh = abs(sim_imm.VelocityVehicle_ms(plot_selector))*3.6;
    abs_v_vehicle_kmh = abs(v_vehicle_kmh);
    [~, v_max_idx] = max(abs_v_vehicle_kmh);
    v_max = sign(v_vehicle_kmh(v_max_idx))*abs_v_vehicle_kmh(v_max_idx);
    v_max = 1.05*v_max;
    v_y_limits = sort([0 v_max],'asc');
    ylim(v_y_limits)
    vfill(mu_geometry_border_times{1},'b','facealpha',.2,'edgecolor','none');
    vfill(mu_geometry_border_times{2},'g','facealpha',.2,'edgecolor','none');
    vfill(mu_geometry_border_times{3},'r','facealpha',.2,'edgecolor','none'); 
    xlabel('t [s]')
    ylabel('v_{ground} [km/h]') 
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','none')
    axis tight    
    ylim(v_y_limits)    
        
    sgtitle(sprintf([ ... 
            'Geometry Identification II\n', ... 
            '(Data Set: ',input_data_selector,'; ', ... 
            '    Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
            '      End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
          ]));
    linkaxes([ax1,ax2,ax3],'x');
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/geom_det_',input_data_selector,'_probab.tex'], ... 
                    'dataPath',['plots/paper_plots/data/geom_det_',input_data_selector,'_probab'], ...
                    'relativeDataPath',['figs/data/geom_det_',input_data_selector,'_probab'], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if

end % if

%% Absolute Error in Cross-Track Direction

if(1)    
                       
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Plot ________________________________________________________________
    
    figure_name = ['Absolute Error in Cross-Track Direction I (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);    
    clear h_plot
    
    clear h_plot    
    h_plot = gobjects(0);   
    ax1 = subplot(2,1,1); hold on; grid on;    
    h_plot_uncertainties = gobjects(0);   
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_time,abs(sim_gnss_filter_input_error_data.e_pp),'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    h_plot(end+1) = plot(sim_ekf_error_data.e_time,abs(sim_ekf_error_data.e_pp),'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.e_time,abs(sim_imm_error_data.e_pp),'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.e_time,abs(sim_tgc_error_data.e_pp),'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    %
    h_plot_uncertainties(end+1) = plot(sim_gnss_filter_input_error_data.e_time,sim_gnss_filter_input_error_data.s.sigma_ct,'k--','LineWidth',1.5,'MarkerSize',10,'DisplayName',sprintf('GPS %.2f Sigma-CT (appr. %.0f%%)',sim_gnss_filter_input_error_data.ee.sigma_gain,sim_gnss_filter_input_error_data.ee.confidence*100));
    h_plot_uncertainties(end+1) = plot(sim_ekf_error_data.e_time,sim_ekf_error_data.s.sigma_ct,'m--','LineWidth',1.5,'MarkerSize',10,'DisplayName',sprintf('EKF %.2f Sigma-CT (appr. %.0f%%)',sim_ekf_error_data.ee.sigma_gain,sim_ekf_error_data.ee.confidence*100));
    h_plot_uncertainties(end+1) = plot(sim_imm_error_data.e_time,sim_imm_error_data.s.sigma_ct,'g--','LineWidth',1.5,'MarkerSize',10,'DisplayName',sprintf('IMM %.2f Sigma-CT (appr. %.0f%%)',sim_imm_error_data.ee.sigma_gain,sim_imm_error_data.ee.confidence*100));
    h_plot_uncertainties(end+1) = plot(sim_tgc_error_data.e_time,sim_tgc_error_data.s.sigma_ct,'b--','LineWidth',1.5,'MarkerSize',10,'DisplayName',sprintf('TGC %.2f Sigma-CT (appr. %.0f%%)',sim_tgc_error_data.ee.sigma_gain,sim_tgc_error_data.ee.confidence*100));
    
    title('Error over time');
    h_legend = legend(h_plot);
    set(h_legend,'Location','northeast')
    xlabel('t [s]')
    ylabel('abs(e_{CT}) [m]')
    axis tight
    
    clear h_plot    
    h_plot = gobjects(0);   
    ax2 = subplot(2,1,2); hold on; grid on;
    plot_selector = true(size(sim_imu_filter_input.Time));
    plot_indices = find(plot_selector);
    h_plot(end+1) = plot(sim_imu_filter_input.Time(plot_selector)-sim_imm.Time(plot_indices(1)),abs(sim_imm.VelocityVehicle_ms(plot_selector))*3.6,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','v_{ground}');
    v_vehicle_kmh = abs(sim_imm.VelocityVehicle_ms(plot_selector))*3.6;
    abs_v_vehicle_kmh = abs(v_vehicle_kmh);
    [~, v_max_idx] = max(abs_v_vehicle_kmh);
    v_max = sign(v_vehicle_kmh(v_max_idx))*abs_v_vehicle_kmh(v_max_idx);
    v_max = 1.05*v_max;
    v_y_limits = sort([0 v_max],'asc');
    axis tight
    ylim(v_y_limits)
    xlabel('t [s]')
    ylabel('v_{ground} [km/h]') 
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','none')
            
        
    title_str = sprintf(['Absolute Perpendicular deviation to reference-map\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(1)),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(end)),')' ... 
                       ]);
    sgtitle(title_str);
    
    linkaxes([ax1,ax2],'x');
    
    % Output __________________________________________________________________

    fprintf('\n### Mean Errors ###\n\n')    
    fprintf('GPS: %.2f\n',mean(abs(sim_gnss_filter_input_error_data.e_pp(~isnan(sim_gnss_filter_input_error_data.e_pp)))));
    fprintf('EKF: %.2f\n',mean(abs(sim_ekf_error_data.e_pp(~isnan(sim_ekf_error_data.e_pp)))));
    fprintf('IMM: %.2f\n',mean(abs(sim_imm_error_data.e_pp(~isnan(sim_imm_error_data.e_pp)))));
    fprintf('TGC: %.2f\n',mean(abs(sim_tgc_error_data.e_pp(~isnan(sim_tgc_error_data.e_pp)))));
    
    fprintf('\n### Max Errors ###\n\n') 
    fprintf('GPS: %.2f\n',max(abs(sim_gnss_filter_input_error_data.e_pp(~isnan(sim_gnss_filter_input_error_data.e_pp)))));
    fprintf('EKF: %.2f\n',max(abs(sim_ekf_error_data.e_pp(~isnan(sim_ekf_error_data.e_pp)))));
    fprintf('IMM: %.2f\n',max(abs(sim_imm_error_data.e_pp(~isnan(sim_imm_error_data.e_pp)))));
    fprintf('TGC: %.2f\n',max(abs(sim_tgc_error_data.e_pp(~isnan(sim_tgc_error_data.e_pp)))));
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/pos_abs_acc_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/pos_abs_acc_',input_data_selector], ...
                    'relativeDataPath',['figs/data/pos_abs_acc_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
    % _____________________________________________________________________
    
    figure_name = ['Absolute Error in Cross-Track Direction II (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);    
    clear h_plot
    
    clear h_plot    
    h_plot = gobjects(0); hold on; grid on; 
    h_plot_uncertainty = gobjects(0);    
    
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_pp_cdf.MaxError_m,sim_gnss_filter_input_error_data.e_pp_cdf.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
    h_plot(end+1) = plot(sim_ekf_error_data.e_pp_cdf.MaxError_m,sim_ekf_error_data.e_pp_cdf.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.e_pp_cdf.MaxError_m,sim_imm_error_data.e_pp_cdf.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.e_pp_cdf.MaxError_m,sim_tgc_error_data.e_pp_cdf.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    h_plot(end+1) = plot(optim_map_error_data.e_pp_cdf.MaxError_m,optim_map_error_data.e_pp_cdf.Availability,'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');
    %
    h_plot_uncertainty(end+1) = plot(sim_gnss_filter_input_error_data.cdf_ct.MaxError_m,sim_gnss_filter_input_error_data.cdf_ct.Availability,'k--','LineWidth',1.5,'MarkerSize',10,'DisplayName','sigma_{GPS}');
    h_plot_uncertainty(end+1) = plot(sim_ekf_error_data.cdf_ct.MaxError_m,sim_ekf_error_data.cdf_ct.Availability,'m--','LineWidth',1.5,'MarkerSize',10,'DisplayName','sigma_{EKF}');
    h_plot_uncertainty(end+1) = plot(sim_imm_error_data.cdf_ct.MaxError_m,sim_imm_error_data.cdf_ct.Availability,'g--','LineWidth',1.5,'MarkerSize',10,'DisplayName','sigma_{IMM}');
    h_plot_uncertainty(end+1) = plot(sim_tgc_error_data.cdf_ct.MaxError_m,sim_tgc_error_data.cdf_ct.Availability,'b--','LineWidth',1.5,'MarkerSize',10,'DisplayName','sigma_{TGC}');
    
    title_str = sprintf([ ... 
                          'CDF\n', ...                           
                          'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h; ', ... 
                          'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'; ', ... 
                          'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)),')' ... 
                       ]);                   
    title(title_str);
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel('error [m]')
    ylabel('CDF')
    axis tight
       
    title_str = sprintf(['Absolute Perpendicular deviation to reference-map\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(1)),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(end)),')' ... 
                       ]);
    sgtitle(title_str);
    
    % Output ______________________________________________________________

    [~,gps_50_idx] = min(abs(sim_gnss_filter_input_error_data.e_pp_cdf.Availability-0.5));
    [~,ekf_50_idx] = min(abs(sim_ekf_error_data.e_pp_cdf.Availability-0.5));
    [~,imm_50_idx] = min(abs(sim_imm_error_data.e_pp_cdf.Availability-0.5));
    [~,tgc_50_idx] = min(abs(sim_tgc_error_data.e_pp_cdf.Availability-0.5));
    
    fprintf('\n### Median Availability Error ###\n\n')    
    fprintf('GPS: %.2f\n',sim_gnss_filter_input_error_data.e_pp_cdf.MaxError_m(gps_50_idx));
    fprintf('EKF: %.2f\n',sim_ekf_error_data.e_pp_cdf.MaxError_m(ekf_50_idx));
    fprintf('IMM: %.2f\n',sim_imm_error_data.e_pp_cdf.MaxError_m(imm_50_idx));
    fprintf('TGC: %.2f\n',sim_tgc_error_data.e_pp_cdf.MaxError_m(tgc_50_idx));
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/pos_stat_acc_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/pos_stat_acc_',input_data_selector], ...
                    'relativeDataPath',['figs/data/pos_stat_acc_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
end % if

%% Availability

if(1)

    % Settings ____________________________________________________________

    positioning_accuracy = 3;

    % Init ________________________________________________________________
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Calculate availabilities at specific localization uncertainty _______

    % GPS 
    cdf_index = find(sim_gnss_filter_input_error_data.cdf_ct.MaxError_m<=positioning_accuracy,1,'last');
    gps_p_availability = sim_gnss_filter_input_error_data.cdf_ct.Availability(cdf_index);

    % EKF
    cdf_index = find(sim_ekf_error_data.cdf_ct.MaxError_m<=positioning_accuracy,1,'last');
    ekf_p_availability = sim_ekf_error_data.cdf_ct.Availability(cdf_index);

    % IMM
    cdf_index = find(sim_imm_error_data.cdf_ct.MaxError_m<=positioning_accuracy,1,'last');
    imm_p_availability = sim_imm_error_data.cdf_ct.Availability(cdf_index);

    % TGC
    cdf_index = find(sim_tgc_error_data.cdf_ct.MaxError_m<=positioning_accuracy,1,'last');
    tgc_p_availability = sim_tgc_error_data.cdf_ct.Availability(cdf_index);

    % Output __________________________________________________________________

    fprintf('\n### Availabilites at abs(e_pos) < %.2f and C = %.2f ###\n\n',positioning_accuracy,error_ellipse_confidence)
    fprintf('GPS availability used in simulation: %.2f\n\n',sim_gnss_filter_input_error_data.availability);
    fprintf('GPS: %.2f\n',gps_p_availability);
    fprintf('EKF: %.2f\n',ekf_p_availability);
    fprintf('IMM: %.2f\n',imm_p_availability);
    fprintf('TGC: %.2f\n',tgc_p_availability);

end % if

%% Cross-Track Error per Availability

if(1)

    % Settings ____________________________________________________________

    positioning_availability = 0.99;

    % Init ________________________________________________________________
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Calculate availabilities at specific localization uncertainty _______

    % GPS 
    cdf_index = find(sim_gnss_filter_input_error_data.cdf_ct.Availability >= positioning_availability,1,'first');
    gps_p_ct_error = sim_gnss_filter_input_error_data.cdf_ct.MaxError_m(cdf_index);

    % EKF
    cdf_index = find(sim_ekf_error_data.cdf_ct.Availability >= positioning_availability,1,'first');
    ekf_p_ct_error = sim_ekf_error_data.cdf_ct.MaxError_m(cdf_index);

    % IMM
    cdf_index = find(sim_imm_error_data.cdf_ct.Availability >= positioning_availability,1,'first');
    imm_p_ct_error = sim_imm_error_data.cdf_ct.MaxError_m(cdf_index);

    % TGC
    cdf_index = find(sim_tgc_error_data.cdf_ct.Availability >= positioning_availability,1,'first');
    tgc_p_ct_error = sim_tgc_error_data.cdf_ct.MaxError_m(cdf_index);

    % Output __________________________________________________________________

    fprintf('\n### CT-Error at V >= %.2f and C = %.2f ###\n\n',positioning_availability,error_ellipse_confidence)
    fprintf('GPS availability used in simulation: %.2f\n\n',sim_gnss_filter_input_error_data.availability);
    fprintf('GPS: %.2f\n',gps_p_ct_error);
    fprintf('EKF: %.2f\n',ekf_p_ct_error);
    fprintf('IMM: %.2f\n',imm_p_ct_error);
    fprintf('TGC: %.2f\n',tgc_p_ct_error);

end % if

%% Uncertainty

if(1)    
                       
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Plot ________________________________________________________________
    
    figure_name = ['Uncertainties (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);    
    clear h_plot
    
    title_str = sprintf([ ... 
                          'Localization Uncertainty\n', ... 
                          '(Data Set: ',input_data_selector,'; ', ... 
                          '   Start: ',datestr(sim_imm.UtcTime(1)),'; ', ... 
                          '     End: ',datestr(sim_imm.UtcTime(end)),')\n' ... 
                          '(conf.: ',num2str(error_ellipse_confidence),'%%; ', ... 
                          'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h; ', ... 
                          'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'; ', ... 
                          'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)),')' ... 
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
    axis tight
%     xlim([0 10])
%     ylim([0 1.2])
    
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
    axis tight
%     xlim([0 10])
%     ylim([0 1.2])
    
    clear h_plot
    h_plot = gobjects(0); 
    h_plot_absolute = gobjects(0);
    ax3 = subplot(1,3,3); hold on; grid on; 
    %title(['CT (conf.: ',num2str(error_ellipse_confidence),', abs(v) >',num2str(v_cdf_min*3.6),'km/h)']);
    title('CT')
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.cdf_ct.MaxError_m,sim_gnss_filter_input_error_data.cdf_ct.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GPS');
    h_plot(end+1) = plot(sim_ekf_error_data.cdf_ct.MaxError_m,sim_ekf_error_data.cdf_ct.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.cdf_ct.MaxError_m,sim_imm_error_data.cdf_ct.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.cdf_ct.MaxError_m,sim_tgc_error_data.cdf_ct.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    %
    if ~strcmp(input_data_selector,'NT')
        h_plot_absolute(end+1) = plot(sim_gnss_filter_input_error_data.e_pp_cdf.MaxError_m,sim_gnss_filter_input_error_data.e_pp_cdf.Availability,'k--','LineWidth',1.5,'MarkerSize',10,'DisplayName','e_{GNSS}');
        h_plot_absolute(end+1) = plot(sim_ekf_error_data.e_pp_cdf.MaxError_m,sim_ekf_error_data.e_pp_cdf.Availability,'m--','LineWidth',1.5,'MarkerSize',10,'DisplayName','e_{EKF}');
        h_plot_absolute(end+1) = plot(sim_imm_error_data.e_pp_cdf.MaxError_m,sim_imm_error_data.e_pp_cdf.Availability,'g--','LineWidth',1.5,'MarkerSize',10,'DisplayName','e_{IMM}');
        h_plot_absolute(end+1) = plot(sim_tgc_error_data.e_pp_cdf.MaxError_m,sim_tgc_error_data.e_pp_cdf.Availability,'b--','LineWidth',1.5,'MarkerSize',10,'DisplayName','e_{TGC}');
        % h_plot_absolute(end+1) = plot(optim_map_error_data.e_pp_cdf.MaxError_m,optim_map_error_data.e_pp_cdf.Availability,'c--','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');
    end % if
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel('error [m]')
    ylabel('availability') 
    axis tight
%     xlim([0 10])
%     ylim([0 1.2])


    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/uncertainties_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/uncertainties_',input_data_selector], ...
                    'relativeDataPath',['figs/data/uncertainties_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
end % if

%% Uncertainty Position Example

if(1)

    plot_error_ellipses = 1; % enable or disable error-ellipse plots    
    plot_optim_data = 1; % plot optimized railway-map
    plot_ref_map_data = 1; % plot reference map
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            plot_start_time = 1995; 
            plot_end_time = 2045;
        case {'BS'}
            plot_start_time = 750; 
            plot_end_time = 900;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_indices = find(plot_selector);    
            
    x_limits = [min(sim_imm.Longitude_deg(plot_selector)) max(sim_imm.Longitude_deg(plot_selector))];
    y_limits = [min(sim_imm.Latitude_deg(plot_selector)) max(sim_imm.Latitude_deg(plot_selector))];
                    
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
    
    title_str = sprintf(['Positions\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
                       ]);
    sgtitle(title_str);
    
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
    
    if exist('x_limits','var')
        %axis equal
        xlim(x_limits)
        ylim(y_limits)
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
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/err_ellip_exam_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/err_ellip_exam_',input_data_selector], ...
                    'relativeDataPath',['figs/data/err_ellip_exam_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/err_ellip_exam_standalone_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/err_ellip_exam_standalone_',input_data_selector], ...
                    'relativeDataPath',['figs/data/err_ellip_exam_standalone_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',true, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
end % if

%% Uncertainty Growth

if(1) 

    fprintf('\n### Mean uncertainty growth with C = %.2f and v_min >= %.2f m/s ###\n\n',error_ellipse_confidence,v_cdf_min)
    fprintf('GPS availability used in simulation: %.2f\n\n',sim_gnss_filter_input_error_data.availability);
    
    fprintf('Mean:\n');
    fprintf('\tEKF: %.2f m/s\n',mean(sim_ekf_error_data.sigma_growth));
    fprintf('\tIMM: %.2f m/s\n',mean(sim_imm_error_data.sigma_growth));
    fprintf('\tTGC: %.2f m/s\n',mean(sim_tgc_error_data.sigma_growth));
    
    fprintf('\nMax:\n');
    fprintf('\tEKF: %.2f m/s\n',max(sim_ekf_error_data.sigma_growth));
    fprintf('\tIMM: %.2f m/s\n',max(sim_imm_error_data.sigma_growth));
    fprintf('\tTGC: %.2f m/s\n',max(sim_tgc_error_data.sigma_growth));

    % figure;histogram(sim_ekf_error_data.sigma_growth,-10:0.1:10);
    % figure;histogram(sim_imm_error_data.sigma_growth,-10:0.1:10);
    % figure;histogram(sim_tgc_error_data.sigma_growth,-10:0.1:10);
    % 
    % [mean(sim_ekf_error_data.sigma_growth), median(sim_ekf_error_data.sigma_growth)]
    % [mean(sim_imm_error_data.sigma_growth), median(sim_imm_error_data.sigma_growth)]
    % [mean(sim_tgc_error_data.sigma_growth), median(sim_tgc_error_data.sigma_growth)]

end % if

%% Initial Map Estimate

if(1)
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            %plot_start_time = 1494; 
            %plot_end_time = 1510;
            plot_start_time = 3655; % 3615; %1494; 
            plot_end_time = 3695; % 3700; % 1510;
        case {'BS'}
            plot_start_time = 750; 
            plot_end_time = 900;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_indices = find(plot_selector);
    
    sim_map_track_element_selector = false(size(sim_map.track_maps.ID));
    for i = 1:length(sim_map.track_maps.ID)
        track_id_i = sim_map.track_maps.ID(i);
        optimization_data_selector_idx = find(track_id_i==optimization_imm_data_selector(:,1));
                
        if any(optimization_imm_data_selector(optimization_data_selector_idx,2:end) & plot_selector(:)')
            sim_map_track_element_selector(i) = true;
        end % if
    end % for i
    
    tgc_map_track_element_selector = false(size(sim_tgc_map.track_id));
    for i = 1:length(sim_tgc_map.track_id)
        track_id_i = sim_tgc_map.track_id(i);
        optimization_data_selector_idx = find(track_id_i==optimization_imm_data_selector(:,1));
                
        if any(optimization_imm_data_selector(optimization_data_selector_idx,2:end) & plot_selector(:)')
            tgc_map_track_element_selector(i) = true;
        end % if
    end % for i
    
    straight_selector_temp = ismember(sim_map_straights.track_start_points.ID,sim_tgc_map.track_id(tgc_map_track_element_selector)); 
    sim_map_straights_temp.topology = sim_map_straights.topology(straight_selector_temp,straight_selector_temp);
    sim_map_straights_temp.track_start_points = sim_map_straights.track_start_points(straight_selector_temp,:);
    sim_map_straights_temp.track_maps = sim_map_straights.track_maps(straight_selector_temp,:);    
    [~,~,sim_map_straights_temp_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_straights_temp,sim_map_density);
    if ~isempty(sim_map_straights_temp_abs_position_utm)
        sim_map_straights_temp_abs_position_utm = p_0_utm(1:2) + sim_map_straights_temp_abs_position_utm;
        [sim_map_straights_temp_abs_position_lat,sim_map_straights_temp_abs_position_lon] = ... 
            utm2ll(sim_map_straights_temp_abs_position_utm(1,:),sim_map_straights_temp_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    else
        sim_map_straights_temp_abs_position_lat = nan;
        sim_map_straights_temp_abs_position_lon = nan;
        
    end % if
    
    arc_selector_temp = ismember(sim_map_arcs.track_start_points.ID,sim_tgc_map.track_id(tgc_map_track_element_selector)); 
    sim_map_arc_temp.topology = sim_map_arcs.topology(arc_selector_temp,arc_selector_temp);
    sim_map_arc_temp.track_start_points = sim_map_arcs.track_start_points(arc_selector_temp,:);
    sim_map_arc_temp.track_maps = sim_map_arcs.track_maps(arc_selector_temp,:);    
    [~,~,sim_map_arc_temp_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_arc_temp,sim_map_density);
    sim_map_arc_temp_abs_position_utm = p_0_utm(1:2) + sim_map_arc_temp_abs_position_utm;
    [sim_map_arcs_temp_abs_position_lat,sim_map_arcs_temp_abs_position_lon] = ... 
        utm2ll(sim_map_arc_temp_abs_position_utm(1,:),sim_map_arc_temp_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
    optim_track_ids = sim_map.track_start_points.ID(sim_map_track_element_selector); 
    [proper_initial_map,~,~] = createProperRailwayMap(sim_map,optimization_data_selector,optim_track_ids);
    
    unkown_selector_temp = (proper_initial_map.track_maps.track_element ~= 1) & ... 
                           (proper_initial_map.track_maps.track_element ~= 3);
    sim_map_unkown_temp.topology = proper_initial_map.topology(unkown_selector_temp,unkown_selector_temp);
    sim_map_unkown_temp.track_start_points = proper_initial_map.track_start_points(unkown_selector_temp,:);
    sim_map_unkown_temp.track_maps = proper_initial_map.track_maps(unkown_selector_temp,:);  
    [~,~,sim_map_unkown_temp_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_unkown_temp,sim_map_density);
    sim_map_unkown_temp_abs_position_utm = p_0_utm(1:2) + sim_map_unkown_temp_abs_position_utm;
    [sim_map_unkown_temp_abs_position_lat,sim_map_unkown_temp_abs_position_lon] = ... 
        utm2ll(sim_map_unkown_temp_abs_position_utm(1,:),sim_map_unkown_temp_abs_position_utm(2,:),p_0_utm(3),'wgs84');
        
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Plot ________________________________________________________________
    
    figure_name = ['Initial Map Estimate (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);      
    clear h_plot    
        
    clear h_plot    
    h_plot = gobjects(0); hold on; grid on;      
    if plot_ref_map
        h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector),sim_gnss_filter_input.Latitude_deg(plot_gnss_selector),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
        axis tight
        x_limits = h_plot(end).Parent.XLim;
        y_limits = h_plot(end).Parent.YLim;
        h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'k-','LineWidth',5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot(end+1) = plot(sim_map_straights_temp_abs_position_lon,sim_map_straights_temp_abs_position_lat,'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
    h_plot(end+1) = plot(sim_map_arcs_temp_abs_position_lon,sim_map_arcs_temp_abs_position_lat,'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');
    h_plot(end+1) = plot(sim_map_unkown_temp_abs_position_lon,sim_map_unkown_temp_abs_position_lat,'r-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (unkown)');
    
%     h_plot(end+1) = plot(optim_map_abs_position_longitude,optim_map_abs_position_latitude,'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');
    
    if exist('x_limits','var')
        %axis equal
        xlim(x_limits)
        ylim(y_limits)
    end % if    
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','southwest')
    xlabel(['Longitude in deg'])
    ylabel(['Latitude in deg'])
%     xlabel(['L [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
%     ylabel(['y [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
    
    % Satellite plot
    % See: https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
    if exist('plot_google_map')
        plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
    end % if
    
    sgtitle(sprintf([ ... 
            'Initial Track-Map\n', ... 
            '(Data Set: ',input_data_selector,'; ', ... 
            '   Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
            '     End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
          ]));
      
    % Output ______________________________________________________________
        
    map_print_format = '%-10i %-10i %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f\n';
    title_data = {'ID','Element','L','r_start','r_end','x_0','y_0','phi_0'};
    print_data = [ ...
                   sim_map.track_maps(sim_map_track_element_selector,:).ID, ...
                   sim_map.track_maps(sim_map_track_element_selector,:).track_element, ... 
                   sim_map.track_maps(sim_map_track_element_selector,:).length, ... 
                   sim_map.track_maps(sim_map_track_element_selector,:).r_start, ... 
                   sim_map.track_maps(sim_map_track_element_selector,:).r_end, ... 
                   sim_map.track_start_points(sim_map_track_element_selector,:).x_0, ... 
                   sim_map.track_start_points(sim_map_track_element_selector,:).y_0, ... 
                   90-sim_map.track_start_points(sim_map_track_element_selector,:).phi_0 ... 
                 ];                  
    
    fprintf('\n### Selected Track-Map ###\n')
    fprintf('%-10s ',title_data{:})
    fprintf('\n')
    fprintf(map_print_format,print_data')
    
    fprintf('\n### UTM-Offset ###\n')
    fprintf('%10.2f\n',p_0_utm)
      
	% matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/initial_map_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/initial_map_',input_data_selector], ...
                    'relativeDataPath',['figs/data/initial_map_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
      
end % if

%% Concatenated Map Estimate

if(1)
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            %plot_start_time = 1494; 
            %plot_end_time = 1510;
            plot_start_time = 3655; % 3615; %1494; 
            plot_end_time = 3695; % 3700; % 1510;
        case {'BS'}
            plot_start_time = 750; 
            plot_end_time = 900;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_indices = find(plot_selector);
    
    sim_map_track_element_selector = false(size(sim_map.track_maps.ID));
    for i = 1:length(sim_map.track_maps.ID)
        track_id_i = sim_map.track_maps.ID(i);
        optimization_data_selector_idx = find(track_id_i==optimization_imm_data_selector(:,1));
                
        if any(optimization_imm_data_selector(optimization_data_selector_idx,2:end) & plot_selector(:)')
            sim_map_track_element_selector(i) = true;
        end % if
    end % for i
    
    optim_track_ids = sim_map.track_start_points.ID(sim_map_track_element_selector); 
    [proper_initial_map,~,~] = createProperRailwayMap(sim_map,optimization_data_selector,optim_track_ids);    
    proper_initial_map.track_start_points{2:end,2:4} = nan(size(proper_initial_map.track_start_points,1)-1,3);
    proper_initial_map.track_start_points = ... 
        calcTableTrackStartPoints([],proper_initial_map);
    proper_initial_map.track_start_points.ID = optim_track_ids;
    proper_initial_map.track_maps.ID = optim_track_ids;
  
    straight_tm_selector = (proper_initial_map.track_maps.track_element == 1);
    straight_sp_selector = ismember(proper_initial_map.track_start_points.ID,unique(proper_initial_map.track_maps(straight_tm_selector,:).ID));
    proper_initial_map_straights.topology = proper_initial_map.topology(straight_sp_selector,straight_sp_selector);
    proper_initial_map_straights.track_start_points = proper_initial_map.track_start_points(straight_sp_selector,:);
    proper_initial_map_straights.track_maps = proper_initial_map.track_maps(straight_tm_selector,:);
    [~,~,proper_initial_map_straights_abs_position_utm,~,~,~,~,~] = calcMapProperties(proper_initial_map_straights,sim_map_density);
    [proper_initial_map_straights_abs_position_lla_latitude,proper_initial_map_straights_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+proper_initial_map_straights_abs_position_utm(1,:),p_0_utm(2)+proper_initial_map_straights_abs_position_utm(2,:),p_0_utm(3),'wgs84');

    arc_tm_selector = (proper_initial_map.track_maps.track_element == 3);
    arc_sp_selector = ismember(proper_initial_map.track_start_points.ID,unique(proper_initial_map.track_maps(arc_tm_selector,:).ID));
    proper_initial_map_arcs.topology = proper_initial_map.topology(arc_sp_selector,arc_sp_selector);
    proper_initial_map_arcs.track_start_points = proper_initial_map.track_start_points(arc_sp_selector,:);
    proper_initial_map_arcs.track_maps = proper_initial_map.track_maps(arc_tm_selector,:);
    [~,~,proper_initial_map_arcs_abs_position_utm,~,~,~,~,~] = calcMapProperties(proper_initial_map_arcs,sim_map_density);
    [proper_initial_map_arcs_abs_position_lla_latitude,proper_initial_map_arcs_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+proper_initial_map_arcs_abs_position_utm(1,:),p_0_utm(2)+proper_initial_map_arcs_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
    unkown_tm_selector = ( ~(proper_initial_map.track_maps.track_element == 3) & ~(proper_initial_map.track_maps.track_element == 1) );
    unkown_sp_selector = ismember(proper_initial_map.track_start_points.ID,unique(proper_initial_map.track_maps(unkown_tm_selector,:).ID));
    proper_initial_map_unkown.topology = proper_initial_map.topology(unkown_sp_selector,unkown_sp_selector);
    proper_initial_map_unkown.track_start_points = proper_initial_map.track_start_points(unkown_sp_selector,:);
    proper_initial_map_unkown.track_maps = proper_initial_map.track_maps(unkown_tm_selector,:);
    [~,~,proper_initial_map_unkown_abs_position_utm,~,~,~,~,~] = calcMapProperties(proper_initial_map_unkown,sim_map_density);
    [proper_initial_map_unkown_abs_position_lla_latitude,proper_initial_map_unkown_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+proper_initial_map_unkown_abs_position_utm(1,:),p_0_utm(2)+proper_initial_map_unkown_abs_position_utm(2,:),p_0_utm(3),'wgs84');
        
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Plot ________________________________________________________________
    
    figure_name = ['Concatenated Map Estimate (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);      
    clear h_plot    
        
    clear h_plot    
    h_plot = gobjects(0); hold on; grid on;      
    if plot_ref_map
        h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector),sim_gnss_filter_input.Latitude_deg(plot_gnss_selector),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
        axis tight
        x_limits = h_plot(end).Parent.XLim;
        y_limits = h_plot(end).Parent.YLim;
        h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'k-','LineWidth',5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot(end+1) = plot(proper_initial_map_straights_abs_position_lla_longitude,proper_initial_map_straights_abs_position_lla_latitude,'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','Concatenated Tracks (straights)');
    h_plot(end+1) = plot(proper_initial_map_arcs_abs_position_lla_longitude,proper_initial_map_arcs_abs_position_lla_latitude,'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','Concatenated Tracks (arcs)');
    h_plot(end+1) = plot(proper_initial_map_unkown_abs_position_lla_longitude,proper_initial_map_unkown_abs_position_lla_latitude,'r-','LineWidth',2.5,'MarkerSize',10,'DisplayName','Concatenated Tracks (unkown)');
    
    %h_plot(end+1) = plot(proper_initial_map_abs_position_longitude,proper_initial_map_abs_position_latitude,'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');
    
    % Annotation
    nan_idx = find(isnan(proper_initial_map_straights_abs_position_lla_longitude));
    mean_start_interval_idx = [1,nan_idx(1:end-1)+1];
    mean_end_interval_idx = nan_idx-1;
    mean_intervals = [mean_start_interval_idx',mean_end_interval_idx'];
    mean_straight_pos_lon = nan(size(mean_intervals,1),1);
    mean_straight_pos_lat = nan(size(mean_intervals,1),1);
    for i = 1:size(mean_intervals,1)
        mean_interval_idx = round(mean([mean_intervals(i,1):mean_intervals(i,2)]));
        mean_straight_pos_lon(i) = (proper_initial_map_straights_abs_position_lla_longitude(mean_interval_idx));
        mean_straight_pos_lat(i) = (proper_initial_map_straights_abs_position_lla_latitude(mean_interval_idx));
    end % for i   
    text(mean_straight_pos_lon,mean_straight_pos_lat,num2str(proper_initial_map_straights.track_maps.ID))
    
    nan_idx = find(isnan(proper_initial_map_arcs_abs_position_lla_longitude));
    mean_start_interval_idx = [1,nan_idx(1:end-1)+1];
    mean_end_interval_idx = nan_idx-1;
    mean_intervals = [mean_start_interval_idx',mean_end_interval_idx'];
    mean_arc_pos_lon = nan(size(mean_intervals,1),1);
    mean_arc_pos_lat = nan(size(mean_intervals,1),1);
    for i = 1:size(mean_intervals,1)
        mean_interval_idx = round(mean([mean_intervals(i,1):mean_intervals(i,2)]));
        mean_arc_pos_lon(i) = (proper_initial_map_arcs_abs_position_lla_longitude(mean_interval_idx));
        mean_arc_pos_lat(i) = (proper_initial_map_arcs_abs_position_lla_latitude(mean_interval_idx));
    end % for i   
    text(mean_arc_pos_lon,mean_arc_pos_lat,num2str(proper_initial_map_arcs.track_maps.ID))
    
    nan_idx = find(isnan(proper_initial_map_unkown_abs_position_lla_longitude));
    mean_start_interval_idx = [1,nan_idx(1:end-1)+1];
    mean_end_interval_idx = nan_idx-1;
    mean_intervals = [mean_start_interval_idx',mean_end_interval_idx'];
    mean_unkown_pos_lon = nan(size(mean_intervals,1),1);
    mean_unkown_pos_lat = nan(size(mean_intervals,1),1);
    for i = 1:size(mean_intervals,1)
        mean_interval_idx = round(mean([mean_intervals(i,1):mean_intervals(i,2)]));
        mean_unkown_pos_lon(i) = (proper_initial_map_unkown_abs_position_lla_longitude(mean_interval_idx));
        mean_unkown_pos_lat(i) = (proper_initial_map_unkown_abs_position_lla_latitude(mean_interval_idx));
    end % for i   
    text(mean_unkown_pos_lon,mean_unkown_pos_lat,num2str(proper_initial_map_unkown.track_maps.ID))
    
    if exist('x_limits','var')
        %axis equal
        xlim(x_limits)
        ylim(y_limits)
    end % if    
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel(['Longitude in deg'])
    ylabel(['Latitude in deg'])
%     xlabel(['L [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
%     ylabel(['y [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
    
    % Satellite plot
    % See: https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
%     if exist('plot_google_map')
%         plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
%     end % if
    
    sgtitle(sprintf([ ... 
            'Concatenated Track-Map\n', ... 
            '(Data Set: ',input_data_selector,'; ', ... 
            '   Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
            '     End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
          ]));
      
    % Output ______________________________________________________________
    
    map_print_format = '%-10i %-10i %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f %-10.1f\n';
    title_data = {'ID','Element','L','r_start','r_end','x_0','y_0','phi_0'};
    print_data = [ ...
                   proper_initial_map.track_maps.ID, ...
                   proper_initial_map.track_maps.track_element, ... 
                   proper_initial_map.track_maps.length, ... 
                   proper_initial_map.track_maps.r_start, ... 
                   proper_initial_map.track_maps.r_end, ... 
                   proper_initial_map.track_start_points.x_0, ... 
                   proper_initial_map.track_start_points.y_0, ... 
                   90-proper_initial_map.track_start_points.phi_0 ... 
                 ];                  
    
    fprintf('\n### Selected Track-Map ###\n')
    fprintf('%-10s ',title_data{:})
    fprintf('\n')
    fprintf(map_print_format,print_data')
    
    fprintf('\n### UTM-Offset ###\n')
    fprintf('%10.2f\n',p_0_utm)
      
	% matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/conca_map_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/conca_map_',input_data_selector], ...
                    'relativeDataPath',['figs/data/conca_map_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
      
end % if

%% Optimized Map Estimate

if(1)
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            %plot_start_time = 1494; 
            %plot_end_time = 1510;
            plot_start_time = 3655; % 3615; %1494; 
            plot_end_time = 3695; % 3700; % 1510;
        case {'BS'}
            plot_start_time = 750; 
            plot_end_time = 900;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_indices = find(plot_selector);    
    
    straight_tm_selector = (optim_map.track_maps.track_element == 1);
    straight_sp_selector = ismember(optim_map.track_start_points.ID,unique(optim_map.track_maps(straight_tm_selector,:).ID));
    optim_map_straights.topology = optim_map.topology(straight_sp_selector,straight_sp_selector);
    optim_map_straights.track_start_points = optim_map.track_start_points(straight_sp_selector,:);
    optim_map_straights.track_maps = optim_map.track_maps(straight_tm_selector,:);
    [~,~,optim_map_straights_abs_position_utm,~,~,~,~,~] = calcMapProperties(optim_map_straights,sim_map_density);
    [optim_map_straights_abs_position_lla_latitude,optim_map_straights_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+optim_map_straights_abs_position_utm(1,:),p_0_utm(2)+optim_map_straights_abs_position_utm(2,:),p_0_utm(3),'wgs84');

    arc_tm_selector = (optim_map.track_maps.track_element == 3);
    arc_sp_selector = ismember(optim_map.track_start_points.ID,unique(optim_map.track_maps(arc_tm_selector,:).ID));
    optim_map_arcs.topology = optim_map.topology(arc_sp_selector,arc_sp_selector);
    optim_map_arcs.track_start_points = optim_map.track_start_points(arc_sp_selector,:);
    optim_map_arcs.track_maps = optim_map.track_maps(arc_tm_selector,:);
    [~,~,optim_map_arcs_abs_position_utm,~,~,~,~,~] = calcMapProperties(optim_map_arcs,sim_map_density);
    [optim_map_arcs_abs_position_lla_latitude,optim_map_arcs_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+optim_map_arcs_abs_position_utm(1,:),p_0_utm(2)+optim_map_arcs_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
    unkown_tm_selector = ( ~(optim_map.track_maps.track_element == 3) & ~(optim_map.track_maps.track_element == 1) );
    unkown_sp_selector = ismember(optim_map.track_start_points.ID,unique(optim_map.track_maps(unkown_tm_selector,:).ID));
    optim_map_unkown.topology = optim_map.topology(unkown_sp_selector,unkown_sp_selector);
    optim_map_unkown.track_start_points = optim_map.track_start_points(unkown_sp_selector,:);
    optim_map_unkown.track_maps = optim_map.track_maps(unkown_tm_selector,:);
    [~,~,optim_map_unkown_abs_position_utm,~,~,~,~,~] = calcMapProperties(optim_map_unkown,sim_map_density);
    [optim_map_unkown_abs_position_lla_latitude,optim_map_unkown_abs_position_lla_longitude] = ... 
        utm2ll(p_0_utm(1)+optim_map_unkown_abs_position_utm(1,:),p_0_utm(2)+optim_map_unkown_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
       
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareRefMapData

    % Plot ________________________________________________________________
    
    figure_name = ['Opt Map Estimate (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);      
    clear h_plot    
        
    clear h_plot    
    h_plot = gobjects(0); hold on; grid on;      
    if plot_ref_map
        h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector),sim_gnss_filter_input.Latitude_deg(plot_gnss_selector),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
        axis tight
        x_limits = h_plot(end).Parent.XLim;
        y_limits = h_plot(end).Parent.YLim;
        h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'k-','LineWidth',5,'MarkerSize',10,'DisplayName','Track-Map');
    end % if
    h_plot(end+1) = plot(optim_map_straights_abs_position_lla_longitude,optim_map_straights_abs_position_lla_latitude,'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','Optim Tracks (straights)');
    h_plot(end+1) = plot(optim_map_arcs_abs_position_lla_longitude,optim_map_arcs_abs_position_lla_latitude,'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','Optim Tracks (arcs)');
    h_plot(end+1) = plot(optim_map_unkown_abs_position_lla_longitude,optim_map_unkown_abs_position_lla_latitude,'r-','LineWidth',2.5,'MarkerSize',10,'DisplayName','Optim Tracks (unkown)');
    
    %h_plot(end+1) = plot(optim_map_abs_position_longitude,optim_map_abs_position_latitude,'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');
    
    if exist('x_limits','var')
        %axis equal
        xlim(x_limits)
        ylim(y_limits)
    end % if    
    
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel(['Longitude in deg'])
    ylabel(['Latitude in deg'])
%     xlabel(['L [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
%     ylabel(['y [m] (UTM zone: ',num2str(p_0_utm(3)),')'])
    
    % Satellite plot
    % See: https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
    if exist('plot_google_map')
        plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
    end % if
    
    sgtitle(sprintf([ ... 
            'Optimized Track-Map\n', ... 
            '(Data Set: ',input_data_selector,'; ', ... 
            '   Start: ',datestr(sim_imm.UtcTime(plot_indices(1))),'; ', ... 
            '     End: ',datestr(sim_imm.UtcTime(plot_indices(end))),')' ... 
          ]));
      
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/final_map_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/final_map_',input_data_selector], ...
                    'relativeDataPath',['figs/data/final_map_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.15g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
      

end % if

%% Map-Error in Cross-Track Direction

if(1)    
                       
    % Plot-Calculations ___________________________________________________  
    
    initLocalization
    prepareSimOutData
    prepareOptimOutData
    prepareRefMapData    

    % Plot ________________________________________________________________
            
    figure_name = ['Map Error in Cross-Track Direction (',input_data_selector,')'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name);    
    clear h_plot
    
    clear h_plot    
    h_plot = gobjects(0); hold on; grid on; 
    h_plot_uncertainty = gobjects(0);    
    
    h_plot(end+1) = plot(sim_gnss_filter_input_error_data.e_pp_cdf.MaxError_m,sim_gnss_filter_input_error_data.e_pp_cdf.Availability,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
    h_plot(end+1) = plot(sim_ekf_error_data.e_pp_cdf.MaxError_m,sim_ekf_error_data.e_pp_cdf.Availability,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
    h_plot(end+1) = plot(sim_imm_error_data.e_pp_cdf.MaxError_m,sim_imm_error_data.e_pp_cdf.Availability,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
    h_plot(end+1) = plot(sim_tgc_error_data.e_pp_cdf.MaxError_m,sim_tgc_error_data.e_pp_cdf.Availability,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
    h_plot(end+1) = plot(optim_map_error_data.e_pp_cdf.MaxError_m,optim_map_error_data.e_pp_cdf.Availability,'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Optim-Map');
        
    title_str = sprintf([ ... 
                          'CDF\n', ...                           
                          'abs(v) > ',num2str(exclude_stillstand_cdf*v_cdf_min*3.6),'km/h; ', ... 
                          'n_{GPS} = ',num2str(length(sim_gnss_filter_input_error_data.cdf_time)),'; ', ... 
                          'n_{Karte} = ',num2str(length(optim_map_abs_pos(:,~isnan(optim_map_abs_pos(1,:))))),'; ', ... 
                          'n_{Filter} = ',num2str(length(sim_tgc_error_data.cdf_time)),')' ... 
                       ]);                   
    title(title_str);
        
    h_legend = legend(h_plot);
    set(h_legend,'Location','southeast')
    xlabel('error [m]')
    ylabel('CDF')
    axis tight
       
    title_str = sprintf(['Absolute Perpendicular deviation to reference-map\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(1)),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(end)),')' ... 
                       ]);
    sgtitle(title_str);
    
    % Output ______________________________________________________________

    fprintf('\n### Mean Errors ###\n\n')    
    fprintf('GPS: %.2f\n',mean(abs(sim_gnss_filter_input_error_data.e_pp(~isnan(sim_gnss_filter_input_error_data.e_pp)))));
    fprintf('EKF: %.2f\n',mean(abs(sim_ekf_error_data.e_pp(~isnan(sim_ekf_error_data.e_pp)))));
    fprintf('IMM: %.2f\n',mean(abs(sim_imm_error_data.e_pp(~isnan(sim_imm_error_data.e_pp)))));
    fprintf('TGC: %.2f\n',mean(abs(sim_tgc_error_data.e_pp(~isnan(sim_tgc_error_data.e_pp)))));
    
    fprintf('\n### Max Errors ###\n\n') 
    fprintf('GPS: %.2f\n',max(abs(sim_gnss_filter_input_error_data.e_pp(~isnan(sim_gnss_filter_input_error_data.e_pp)))));
    fprintf('EKF: %.2f\n',max(abs(sim_ekf_error_data.e_pp(~isnan(sim_ekf_error_data.e_pp)))));
    fprintf('IMM: %.2f\n',max(abs(sim_imm_error_data.e_pp(~isnan(sim_imm_error_data.e_pp)))));
    fprintf('TGC: %.2f\n',max(abs(sim_tgc_error_data.e_pp(~isnan(sim_tgc_error_data.e_pp)))));
    fprintf('Map: %.2f\n',max(abs(optim_map_error_data.e_pp(~isnan(optim_map_error_data.e_pp)))));
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/mapping_error_',input_data_selector,'.tex'], ... 
                    'dataPath',['plots/paper_plots/data/mapping_error_',input_data_selector], ...
                    'relativeDataPath',['figs/data/mapping_error_',input_data_selector], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
    
end % if

% %% Optim-Map Track-Length Error CDF
% 
% if(1)
%     
%     [~,p_ref_pp,~] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,optim_map_abs_pos);
%     %p_ref_pp = p_ref_pp - p_ref_pp(:,1);
%        
%     d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp,1,2)].^2,1))); 
%     
%    optim_map_temp_abs_position = [];
%    d_p_ref_pp_selector = false(size(d_p_ref_pp));
%     for i = 1:size(optim_map.track_maps,1)
%         track_map_i = optim_map.track_maps(i,:);
%         d_train_i = d_p_ref_pp;
%         d_length = optim_map.track_maps(i,:).length;
%         if i > 1
%             d_before = sum(optim_map.track_maps(1:i-1,:).length);
%    
%             d_train_i = d_p_ref_pp-d_before;
%             
%             d_p_ref_pp_selector( (d_train_i>=0) & (d_train_i <= d_length) ) = true;
%             
%             d_train_i = d_train_i(d_train_i>=0);
%         else
%             d_p_ref_pp_selector(d_p_ref_pp-d_before <= d_length) = true;
%         end % if
%         %d_p_ref_pp_selector = d_p_ref_pp_selector & (d_p_ref_pp-d_before <= d_length);
%         d_train_i = d_train_i(d_train_i<=d_length);
%         
%         if isempty(d_train_i)
%             abs_pos = nan(2,1);
%         else
%             [~,~,~,abs_pos,~,~,~,~,~,~] = calcTrackRouteProperties(d_train_i,track_map_i.ID,track_map_i.ID(1),0,optim_map,1);
%         end % if        
%         
%         optim_map_temp_abs_position = [optim_map_temp_abs_position, abs_pos];
%     end % for i
%     
% %     [~,~,~,optim_map_temp_abs_position,~,~,~,~] = calcTrackProperties(0:1:600,optim_track_maps.ID(1),0,optim_track_maps,optim_map.track_start_points(1,:).phi_0);
% %     optim_map_temp_abs_position = ... 
% %         optim_map_temp_abs_position + [optim_map.track_start_points(1,:).x_0;optim_map.track_start_points(1,:).y_0;];
% %     [~,~,~,optim_map_temp_abs_position,~,~,~,~,~,~] = calcTrackRouteProperties(d_p_ref_pp,optim_track_maps.ID,optim_track_maps.ID(1),0,optim_map,1);
% %     [~,~,optim_map_temp_abs_position,~,~,~,~,~] = calcMapProperties(optim_map,1);
%     
%     e_p_tl = sqrt(sum((p_ref_pp(:,d_p_ref_pp_selector) - optim_map_temp_abs_position).^2,1));
%     plot(d_p_ref_pp(d_p_ref_pp_selector),e_p_tl);
%     
% %     plot(p_ref_pp(1,:),p_ref_pp(2,:)); hold all
% %     plot(optim_map_temp_abs_position(1,:),optim_map_temp_abs_position(2,:));
% %     plot(optim_map_abs_pos(1,:),optim_map_abs_pos(2,:));
% 
% end % if

%% Optim-Map Track-Length Final Error

if(1)
    
    [~,p_ref_pp,p_ref_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,optim_map_abs_pos);
       
    d_p_ref_pp = cumsum(sqrt(sum([zeros(2,1), diff(p_ref_pp,1,2)].^2,1)));    
    p_optim_map_temp = optim_map_abs_pos(:,p_ref_pp_indices);
    d_p_optim_map = cumsum(sqrt(sum([zeros(2,1),diff(p_optim_map_temp,1,2)].^2,1)));
    
    fprintf('\n### Track-Length ###\n\n') 
    fprintf('Reference Map Length: %.2f\n',max(d_p_ref_pp));
    fprintf('Optim Map Length: %.2f\n',max(d_p_optim_map));
    
    fprintf('\n### Diffs (min,max) ###\n\n') 
    fprintf('GPS: %10.1f \n',sim_gnss_filter_input_error_data.e_tl(end));
    fprintf('IMM: %10.1f \n',sim_ekf_error_data.e_tl(end));
    fprintf('EKF: %10.1f \n',sim_imm_error_data.e_tl(end));
    fprintf('TGC: %10.1f \n',sim_tgc_error_data.e_tl(end));
    fprintf('Map: %10.1f \n',optim_map_error_data.e_tl(end));

end % if

%% Perpendicular Error Visualization (GPS to Track-Elements)

if(1)
    
    plot_llh = 0;
    plot_utm = 1;
    plot_pp_error = 1;
    plot_ref_map_data = 0; % plot reference map
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            plot_start_time = 1865; % 3750; % 1800; 
            plot_end_time = 1885; % 3820; % 1930;
        case {'BS'}
            plot_start_time = 10; 
            plot_end_time = 100;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_sim_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_sim_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_sim_indices = find(plot_sim_selector);    
    plot_gnss_indices = find(plot_gnss_selector);
            
    x_llh_limits = [min(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector)) max(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector))];
    y_llh_limits = [min(sim_gnss_filter_input.Latitude_deg(plot_gnss_selector)) max(sim_gnss_filter_input.Latitude_deg(plot_gnss_selector))];
    x_utm_limits = [min(sim_gnss_filter_input.UtmEast_m(plot_gnss_selector)) max(sim_gnss_filter_input.UtmEast_m(plot_gnss_selector))];
    y_utm_limits = [min(sim_gnss_filter_input.UtmNorth_m(plot_gnss_selector)) max(sim_gnss_filter_input.UtmNorth_m(plot_gnss_selector))];     
    
    % Bounded reference map
    ref_map_abs_position_utm = [imm_ref_track_map.UtmEast_m(plot_sim_selector)';imm_ref_track_map.UtmNorth_m(plot_sim_selector)'];
    ref_map_abs_position_ll = [imm_ref_track_map.Latitude_deg(plot_sim_selector)';imm_ref_track_map.Longitude_deg(plot_sim_selector)'];
    
    sim_map_track_element_selector = false(size(sim_map.track_maps.ID));
    for i = 1:length(sim_map.track_maps.ID)
        track_id_i = sim_map.track_maps.ID(i);
        optimization_data_selector_idx = find(track_id_i==optimization_imm_data_selector(:,1));
                
        if any(optimization_imm_data_selector(optimization_data_selector_idx,2:end) & plot_sim_selector(:)')
            sim_map_track_element_selector(i) = true;
        end % if
    end % for i
    
    % Straights PP-Error
    sim_map_straights_temp_abs_position_utm = [];
    sim_map_straights_temp_abs_position_lat = [];
    sim_map_straights_temp_abs_position_lon = [];
    straights_gnss_pos = [];
    p_ref_straights_gnss_pp = [];
    straights_gnss_pos_temp_lat = [];
    straights_gnss_pos_temp_lon = [];
    p_ref_straights_gnss_pp_lat = [];
    p_ref_straights_gnss_pp_lon = [];
    for i = 1:length(sim_map_straights.track_maps.ID)
        id_i = sim_map_straights.track_maps(i,:).ID;        
        optimization_gps_data_selector_i = (optimization_gps_data_selector(:,1) == id_i);
        data_selector_i = optimization_gps_data_selector(optimization_gps_data_selector_i,2:end);
        data_selector_i = data_selector_i(:) & plot_gnss_selector(:);
        if any(data_selector_i)
            sim_map_straights_temp.topology = 0;
            sim_map_straights_temp.track_start_points = sim_map_straights.track_start_points(i,:);
            sim_map_straights_temp.track_maps = sim_map_straights.track_maps(i,:);  
            
            [~,~,sim_map_straights_temp_abs_position_utm_i,~,~,~,~,~] = calcMapProperties(sim_map_straights_temp,sim_map_density);
            sim_map_straights_temp_abs_position_utm_i = p_0_utm(1:2) + sim_map_straights_temp_abs_position_utm_i;
            [sim_map_straights_temp_abs_position_lat_i,sim_map_straights_temp_abs_position_lon_i] = ... 
                utm2ll(sim_map_straights_temp_abs_position_utm_i(1,:),sim_map_straights_temp_abs_position_utm_i(2,:),p_0_utm(3),'wgs84');
            
            sim_map_straights_temp_abs_position_utm = [sim_map_straights_temp_abs_position_utm, sim_map_straights_temp_abs_position_utm_i];
            sim_map_straights_temp_abs_position_lat = [sim_map_straights_temp_abs_position_lat, sim_map_straights_temp_abs_position_lat_i];
            sim_map_straights_temp_abs_position_lon = [sim_map_straights_temp_abs_position_lon, sim_map_straights_temp_abs_position_lon_i];
            
            % GPS PP-Error (...)            
            straights_gnss_pos_i = [ sim_gnss_filter_input.UtmEast_m(data_selector_i)'; sim_gnss_filter_input.UtmNorth_m(data_selector_i)'];
            [~,p_ref_straights_gnss_pp_i,p_ref_straights_gnss_pp_indices_i] = calcPerpendicularErrorToRefPoints(sim_map_straights_temp_abs_position_utm_i,straights_gnss_pos_i);
            straights_gnss_pos_i = straights_gnss_pos_i(:,p_ref_straights_gnss_pp_indices_i);
            [straights_gnss_pos_i_temp_lat,straights_gnss_pos_i_temp_lon] = utm2ll(straights_gnss_pos_i(1,:)',straights_gnss_pos_i(2,:)',p_0_utm(3));
            [p_ref_straights_gnss_pp_lat_i,p_ref_straights_gnss_pp_lon_i] = utm2ll(p_ref_straights_gnss_pp_i(1,:)',p_ref_straights_gnss_pp_i(2,:)',p_0_utm(3));
            
            straights_gnss_pos = [straights_gnss_pos,straights_gnss_pos_i];
            p_ref_straights_gnss_pp = [p_ref_straights_gnss_pp,p_ref_straights_gnss_pp_i];
            straights_gnss_pos_temp_lat = [straights_gnss_pos_temp_lat;straights_gnss_pos_i_temp_lat];
            straights_gnss_pos_temp_lon = [straights_gnss_pos_temp_lon;straights_gnss_pos_i_temp_lon];
            p_ref_straights_gnss_pp_lat = [p_ref_straights_gnss_pp_lat;p_ref_straights_gnss_pp_lat_i];
            p_ref_straights_gnss_pp_lon = [p_ref_straights_gnss_pp_lon;p_ref_straights_gnss_pp_lon_i];
    
        else
            sim_map_straights_temp_abs_position_utm = [sim_map_straights_temp_abs_position_utm, nan(2,1)];
            sim_map_straights_temp_abs_position_lat = [sim_map_straights_temp_abs_position_lat, nan];
            sim_map_straights_temp_abs_position_lon = [sim_map_straights_temp_abs_position_lon, nan];
            
            straights_gnss_pos = [straights_gnss_pos,nan(2,1)];
            p_ref_straights_gnss_pp = [p_ref_straights_gnss_pp,nan(2,1)];
            straights_gnss_pos_temp_lat = [straights_gnss_pos_temp_lat;nan];
            straights_gnss_pos_temp_lon = [straights_gnss_pos_temp_lon;nan];
            p_ref_straights_gnss_pp_lat = [p_ref_straights_gnss_pp_lat;nan];
            p_ref_straights_gnss_pp_lon = [p_ref_straights_gnss_pp_lon;nan];
        end % if
    end % if
    
    % Arcs PP-Error    
    sim_map_arcs_temp_abs_position_utm = [];
    sim_map_arcs_temp_abs_position_lat = [];
    sim_map_arcs_temp_abs_position_lon = [];
    arcs_gnss_pos = [];
    p_ref_arcs_gnss_pp = [];
    arcs_gnss_pos_temp_lat = [];
    arcs_gnss_pos_temp_lon = [];
    p_ref_arcs_gnss_pp_lat = [];
    p_ref_arcs_gnss_pp_lon = [];
    for i = 1:length(sim_map_arcs.track_maps.ID)
        id_i = sim_map_arcs.track_maps(i,:).ID;        
        optimization_gps_data_selector_i = (optimization_gps_data_selector(:,1) == id_i);
        data_selector_i = optimization_gps_data_selector(optimization_gps_data_selector_i,2:end);
        data_selector_i = data_selector_i(:) & plot_gnss_selector(:);
        if any(data_selector_i)
            sim_map_arcs_temp.topology = 0;
            sim_map_arcs_temp.track_start_points = sim_map_arcs.track_start_points(i,:);
            sim_map_arcs_temp.track_maps = sim_map_arcs.track_maps(i,:);  
            
            [~,~,sim_map_arcs_temp_abs_position_utm_i,~,~,~,~,~] = calcMapProperties(sim_map_arcs_temp,sim_map_density);
            sim_map_arcs_temp_abs_position_utm_i = p_0_utm(1:2) + sim_map_arcs_temp_abs_position_utm_i;
            [sim_map_arcs_temp_abs_position_lat_i,sim_map_arcs_temp_abs_position_lon_i] = ... 
                utm2ll(sim_map_arcs_temp_abs_position_utm_i(1,:),sim_map_arcs_temp_abs_position_utm_i(2,:),p_0_utm(3),'wgs84');
            
            sim_map_arcs_temp_abs_position_utm = [sim_map_arcs_temp_abs_position_utm, sim_map_arcs_temp_abs_position_utm_i];
            sim_map_arcs_temp_abs_position_lat = [sim_map_arcs_temp_abs_position_lat, sim_map_arcs_temp_abs_position_lat_i];
            sim_map_arcs_temp_abs_position_lon = [sim_map_arcs_temp_abs_position_lon, sim_map_arcs_temp_abs_position_lon_i];
            
            % GPS PP-Error (...)            
            arcs_gnss_pos_i = [ sim_gnss_filter_input.UtmEast_m(data_selector_i)'; sim_gnss_filter_input.UtmNorth_m(data_selector_i)'];
            [~,p_ref_arcs_gnss_pp_i,p_ref_arcs_gnss_pp_indices_i] = calcPerpendicularErrorToRefPoints(sim_map_arcs_temp_abs_position_utm_i,arcs_gnss_pos_i);
            arcs_gnss_pos_i = arcs_gnss_pos_i(:,p_ref_arcs_gnss_pp_indices_i);
            [arcs_gnss_pos_i_temp_lat,arcs_gnss_pos_i_temp_lon] = utm2ll(arcs_gnss_pos_i(1,:)',arcs_gnss_pos_i(2,:)',p_0_utm(3));
            [p_ref_arcs_gnss_pp_lat_i,p_ref_arcs_gnss_pp_lon_i] = utm2ll(p_ref_arcs_gnss_pp_i(1,:)',p_ref_arcs_gnss_pp_i(2,:)',p_0_utm(3));
            
            arcs_gnss_pos = [arcs_gnss_pos,arcs_gnss_pos_i];
            p_ref_arcs_gnss_pp = [p_ref_arcs_gnss_pp,p_ref_arcs_gnss_pp_i];
            arcs_gnss_pos_temp_lat = [arcs_gnss_pos_temp_lat;arcs_gnss_pos_i_temp_lat];
            arcs_gnss_pos_temp_lon = [arcs_gnss_pos_temp_lon;arcs_gnss_pos_i_temp_lon];
            p_ref_arcs_gnss_pp_lat = [p_ref_arcs_gnss_pp_lat;p_ref_arcs_gnss_pp_lat_i];
            p_ref_arcs_gnss_pp_lon = [p_ref_arcs_gnss_pp_lon;p_ref_arcs_gnss_pp_lon_i];
    
        else
            sim_map_arcs_temp_abs_position_utm = [sim_map_arcs_temp_abs_position_utm, nan(2,1)];
            sim_map_arcs_temp_abs_position_lat = [sim_map_arcs_temp_abs_position_lat, nan];
            sim_map_arcs_temp_abs_position_lon = [sim_map_arcs_temp_abs_position_lon, nan];
            
            arcs_gnss_pos = [arcs_gnss_pos,nan(2,1)];
            p_ref_arcs_gnss_pp = [p_ref_arcs_gnss_pp,nan(2,1)];
            arcs_gnss_pos_temp_lat = [arcs_gnss_pos_temp_lat;nan];
            arcs_gnss_pos_temp_lon = [arcs_gnss_pos_temp_lon;nan];
            p_ref_arcs_gnss_pp_lat = [p_ref_arcs_gnss_pp_lat;nan];
            p_ref_arcs_gnss_pp_lon = [p_ref_arcs_gnss_pp_lon;nan];
        end % if
    end % if
    
    % Unkown elements
    optim_track_ids = sim_map.track_start_points.ID(sim_map_track_element_selector); 
    [proper_initial_map,~,~] = createProperRailwayMap(sim_map,optimization_data_selector,optim_track_ids);
    
    unkown_selector_temp = (proper_initial_map.track_maps.track_element ~= 1) & ... 
                           (proper_initial_map.track_maps.track_element ~= 3);
    sim_map_unkown_temp.topology = proper_initial_map.topology(unkown_selector_temp,unkown_selector_temp);
    sim_map_unkown_temp.track_start_points = proper_initial_map.track_start_points(unkown_selector_temp,:);
    sim_map_unkown_temp.track_maps = proper_initial_map.track_maps(unkown_selector_temp,:);  
    [~,~,sim_map_unkown_temp_abs_position_utm,~,~,~,~,~] = calcMapProperties(sim_map_unkown_temp,sim_map_density);
    sim_map_unkown_temp_abs_position_utm = p_0_utm(1:2) + sim_map_unkown_temp_abs_position_utm;
    [sim_map_unkown_temp_abs_position_lat,sim_map_unkown_temp_abs_position_lon] = ... 
        utm2ll(sim_map_unkown_temp_abs_position_utm(1,:),sim_map_unkown_temp_abs_position_utm(2,:),p_0_utm(3),'wgs84');
    
%     % Optim Map PP-Error    
%     [~,p_ref_optim_map_pp,p_ref_optim_map_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,optim_map_abs_pos+p_0_utm(1:2));
%     optim_map_abs_pos_temp = optim_map_abs_pos(:,p_ref_optim_map_pp_indices)+p_0_utm(1:2);
%     [optim_map_abs_pos_temp_lat,optim_map_abs_pos_temp_lon] = utm2ll(optim_map_abs_pos_temp(1,:)',optim_map_abs_pos_temp(2,:)',p_0_utm(3));
%     [p_ref_optim_map_pp_lat,p_ref_optim_map_pp_lon] = utm2ll(p_ref_optim_map_pp(1,:)',p_ref_optim_map_pp(2,:)',p_0_utm(3)); 
    
%     plot(ref_map_abs_position_utm(1,:),ref_map_abs_position_utm(2,:)); hold on;
%     %plot(sim_map_arcs_temp_abs_position_utm(1,:),sim_map_arcs_temp_abs_position_utm(2,:)); hold on
%     %plot(sim_map_straights_temp_abs_position_utm(1,:),sim_map_straights_temp_abs_position_utm(2,:));
%     for i = 1:length(arcs_gnss_pos)
%         plot([arcs_gnss_pos(1,i),p_ref_arcs_gnss_pp(1,i)],[arcs_gnss_pos(2,i),p_ref_arcs_gnss_pp(2,i)]);
%     end % for i
%     for i = 1:length(straights_gnss_pos)
%         plot([straights_gnss_pos(1,i),p_ref_straights_gnss_pp(1,i)],[straights_gnss_pos(2,i),p_ref_straights_gnss_pp(2,i)]);
%     end % for i
                    
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
    
    figure_name = 'PP-Error Visualization (GPS)';
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold on;
    clear h_plot    
    h_plot = gobjects(0);
    h_d_pp = gobjects(0);
    
    title_str = sprintf(['Positions\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(plot_sim_indices(1))),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(plot_sim_indices(end))),')' ... 
                       ]);
    sgtitle(title_str);
    
    hold on; grid on;    
       
%     h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg,sim_gnss_filter_input.Latitude_deg,'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
%     h_plot(end+1) = plot(sim_imm.Longitude_deg,sim_imm.Latitude_deg,'g-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMM');
%     h_plot(end+1) = plot(sim_ekf.Longitude_deg,sim_ekf.Latitude_deg,'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','EKF');
%     h_plot(end+1) = plot(sim_tgc.Longitude_deg,sim_tgc.Latitude_deg,'b-','LineWidth',1.5,'MarkerSize',10,'DisplayName','TGC');
%     % h_plot(end+1) = plot(sim_tgc_map.Longitude_deg,sim_tgc_map.Latitude_deg,'r-','LineWidth',3,'MarkerSize',10,'DisplayName','TGC Tracks');
%     h_plot(end+1) = plot(sim_tgc_map.Straights_Longitude_deg,sim_tgc_map.Straights_Latitude_deg,'-','Color',[1 0.5 0],'LineWidth',3,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
%     h_plot(end+1) = plot(sim_tgc_map.Arcs_Longitude_deg,sim_tgc_map.Arcs_Latitude_deg,'-','Color',[1 1 1],'LineWidth',3,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');

    if plot_llh
        if plot_pp_error
            for i = 1:length(straights_gnss_pos)
                h_d_pp(end+1) = plot([straights_gnss_pos_temp_lon(i,1),p_ref_straights_gnss_pp_lon(i,1)],[straights_gnss_pos_temp_lat(i,1),p_ref_straights_gnss_pp_lat(i,1)],'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (straights)');
            end % for i
            for i = 1:length(arcs_gnss_pos)
                h_d_pp(end+1) = plot([arcs_gnss_pos_temp_lon(i,1),p_ref_arcs_gnss_pp_lon(i,1)],[arcs_gnss_pos_temp_lat(i,1),p_ref_arcs_gnss_pp_lat(i,1)],'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (arcs)');
            end % for i            
%             for i = 1:length(optim_map_abs_pos_temp)
%                 h_d_pp(end+1) = plot([optim_map_abs_pos_temp_lon(i,1),p_ref_optim_map_pp_lon(i,1)],[optim_map_abs_pos_temp_lat(i,1),p_ref_optim_map_pp_lat(i,1)],'y-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (Optim-Map)');
%             end % for i
        end % if
        
        h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg(plot_gnss_indices),sim_gnss_filter_input.Latitude_deg(plot_gnss_indices),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
        h_plot(end+1) = plot(sim_map_straights_temp_abs_position_lon,sim_map_straights_temp_abs_position_lat,'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
        h_plot(end+1) = plot(sim_map_arcs_temp_abs_position_lon,sim_map_arcs_temp_abs_position_lat,'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');
        h_plot(end+1) = plot(sim_map_unkown_temp_abs_position_lon,sim_map_unkown_temp_abs_position_lat,'r-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (unkown)');
%         h_plot(end+1) = plot(optim_map_abs_position_longitude,optim_map_abs_position_latitude,'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');

        if plot_ref_map_data && plot_ref_map
            h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
        end % if
        
        
        if exist('x_llh_limits','var')
            %axis equal
            xlim(x_llh_limits)
            ylim(y_llh_limits)
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
        
    if plot_utm
        if plot_pp_error
            for i = 1:length(straights_gnss_pos)
                h_d_pp(end+1) = plot([straights_gnss_pos(1,i),p_ref_straights_gnss_pp(1,i)],[straights_gnss_pos(2,i),p_ref_straights_gnss_pp(2,i)],'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (straights)');
            end % for i
            for i = 1:length(arcs_gnss_pos)
                h_d_pp(end+1) = plot([arcs_gnss_pos(1,i),p_ref_arcs_gnss_pp(1,i)],[arcs_gnss_pos(2,i),p_ref_arcs_gnss_pp(2,i)],'m-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (arcs)');
            end % for i
        end % if
        
        h_plot(end+1) = plot(sim_gnss_filter_input.UtmEast_m(plot_gnss_indices),sim_gnss_filter_input.UtmNorth_m(plot_gnss_indices),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
        h_plot(end+1) = plot(sim_map_straights_temp_abs_position_utm(1,:),sim_map_straights_temp_abs_position_utm(2,:),'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
        h_plot(end+1) = plot(sim_map_arcs_temp_abs_position_utm(1,:),sim_map_arcs_temp_abs_position_utm(2,:),'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');
        h_plot(end+1) = plot(sim_map_unkown_temp_abs_position_utm(1,:),sim_map_unkown_temp_abs_position_utm(2,:),'r-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (unkown)');
        
        if plot_ref_map_data && plot_ref_map
            h_plot(end+1) = plot(ref_track_map.UtmEast_m,ref_track_map.UtmNorth_m,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
        end % if 
        
        if exist('x_utm_limits','var')
            axis equal
%             xlim(x_utm_limits)
%             ylim(y_utm_limits)
        end % if  

        h_legend = legend(h_plot);
        set(h_legend,'Location','northeast')
        xlabel('east [m]')
        ylabel('north [m]')
        % axis equal        
    end % if
    
        % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/pp_error_example','.tex'], ... 
                    'dataPath',['plots/paper_plots/data/pp_error_example'], ...
                    'relativeDataPath',['figs/data/pp_error_example'], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
        
end % if

%% Perpendicular Error Visualization (Optim-Map to Ref-Map)

if(0)
    
    plot_llh = 0;
    plot_utm = 1;
    plot_pp_error = 1;
    plot_ref_map_data = 1; % plot reference map
    
    use_time_frame = 1;
    switch input_data_selector
        case {'C'}
            plot_start_time = 0; % 3750; % 1800; 
            plot_end_time = 600; % 3820; % 1930;
        case {'BS'}
            plot_start_time = 10; 
            plot_end_time = 100;
        case {'NT'}
            plot_start_time = 10; 
            plot_end_time = 100;
    end % switch
    
    if use_time_frame
        plot_sim_selector = (sim_imm.Time >= plot_start_time) & (sim_imm.Time <= plot_end_time);
        plot_gnss_selector = (sim_gnss_filter_input.Time >= plot_start_time) & (sim_gnss_filter_input.Time <= plot_end_time);
    else
        plot_sim_selector = true(size(sim_imm.Time));
        plot_gnss_selector = true(size(sim_gnss_filter_input.Time));
    end % if
    plot_sim_indices = find(plot_sim_selector);    
    plot_gnss_indices = find(plot_gnss_selector);
            
    x_llh_limits = [min(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector)) max(sim_gnss_filter_input.Longitude_deg(plot_gnss_selector))];
    y_llh_limits = [min(sim_gnss_filter_input.Latitude_deg(plot_gnss_selector)) max(sim_gnss_filter_input.Latitude_deg(plot_gnss_selector))];
    x_utm_limits = [min(sim_gnss_filter_input.UtmEast_m(plot_gnss_selector)) max(sim_gnss_filter_input.UtmEast_m(plot_gnss_selector))];
    y_utm_limits = [min(sim_gnss_filter_input.UtmNorth_m(plot_gnss_selector)) max(sim_gnss_filter_input.UtmNorth_m(plot_gnss_selector))];     
    
    % Bounded reference map
    ref_map_abs_position_utm = [imm_ref_track_map.UtmEast_m(plot_sim_selector)';imm_ref_track_map.UtmNorth_m(plot_sim_selector)'];
%     ref_map_abs_position_utm = [ref_track_map.UtmEast_m';ref_track_map.UtmNorth_m'];
    ref_map_abs_position_ll = [imm_ref_track_map.Latitude_deg(plot_sim_selector)';imm_ref_track_map.Longitude_deg(plot_sim_selector)'];
            
    % Optim Map PP-Error    
    [~,p_ref_optim_map_pp,p_ref_optim_map_pp_indices] = calcPerpendicularErrorToRefPoints(ref_map_abs_position_utm,optim_map_abs_pos+p_0_utm(1:2));
    optim_map_abs_pos_temp = optim_map_abs_pos(:,p_ref_optim_map_pp_indices)+p_0_utm(1:2);
    [optim_map_abs_pos_temp_lat,optim_map_abs_pos_temp_lon] = utm2ll(optim_map_abs_pos_temp(1,:)',optim_map_abs_pos_temp(2,:)',p_0_utm(3));
    [p_ref_optim_map_pp_lat,p_ref_optim_map_pp_lon] = utm2ll(p_ref_optim_map_pp(1,:)',p_ref_optim_map_pp(2,:)',p_0_utm(3)); 
                        
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
    
    figure_name = 'PP-Error Visualization (Optim-Map)';
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold on;
    clear h_plot    
    h_plot = gobjects(0);
    h_d_pp = gobjects(0);
    
    title_str = sprintf(['Positions\n' ... 
                         '(Data Set: ',input_data_selector,'; ', ... 
                         '   Start: ',datestr(sim_imm.UtcTime(plot_sim_indices(1))),'; ', ... 
                         '     End: ',datestr(sim_imm.UtcTime(plot_sim_indices(end))),')' ... 
                       ]);
    sgtitle(title_str);
    
    hold on; grid on;    
    if plot_llh
        if plot_pp_error         
            for i = 1:length(optim_map_abs_pos_temp)
                h_d_pp(end+1) = plot([optim_map_abs_pos_temp_lon(i,1),p_ref_optim_map_pp_lon(i,1)],[p_ref_optim_map_pp_lat(i,1),p_ref_optim_map_pp_lat(i,1)],'y-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (Optim-Map)');
            end % for i
        end % if
        
%         h_plot(end+1) = plot(sim_gnss_filter_input.Longitude_deg(plot_gnss_indices),sim_gnss_filter_input.Latitude_deg(plot_gnss_indices),'o','LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName','GNSS');
%         h_plot(end+1) = plot(sim_map_straights_temp_abs_position_lon,sim_map_straights_temp_abs_position_lat,'c-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (straights)');
%         h_plot(end+1) = plot(sim_map_arcs_temp_abs_position_lon,sim_map_arcs_temp_abs_position_lat,'m-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (arcs)');
%         h_plot(end+1) = plot(sim_map_unkown_temp_abs_position_lon,sim_map_unkown_temp_abs_position_lat,'r-','LineWidth',2.5,'MarkerSize',10,'DisplayName','TGC Tracks (unkown)');
        h_plot(end+1) = plot(optim_map_abs_position_longitude,optim_map_abs_position_latitude,'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');

        if plot_ref_map_data && plot_ref_map
            h_plot(end+1) = plot(ref_track_map.Longitude_deg,ref_track_map.Latitude_deg,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
        end % if
        
        
        if exist('x_llh_limits','var')
            %axis equal
            xlim(x_llh_limits)
            ylim(y_llh_limits)
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
        
    if plot_utm
        if plot_pp_error
            for i = 1:length(optim_map_abs_pos_temp)
                h_d_pp(end+1) = plot([optim_map_abs_pos_temp(1,i),p_ref_optim_map_pp(1,i)],[optim_map_abs_pos_temp(2,i),p_ref_optim_map_pp(2,i)],'c-','LineWidth',1.5,'MarkerSize',10,'DisplayName','PP-Error (straights)');
            end % for i
        end % if
        
        h_plot(end+1) = plot(optim_map_abs_pos_temp(1,:),optim_map_abs_pos_temp(2,:),'c-','LineWidth',3,'MarkerSize',10,'DisplayName','Optim Map');
        
        if plot_ref_map_data && plot_ref_map
            h_plot(end+1) = plot(ref_track_map.UtmEast_m,ref_track_map.UtmNorth_m,'r-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Track-Map');  
        end % if 
        
        if exist('x_utm_limits','var')
            axis equal
            xlim(x_utm_limits)
            ylim(y_utm_limits)
        end % if  

        h_legend = legend(h_plot);
        set(h_legend,'Location','northeast')
        xlabel('east [m]')
        ylabel('north [m]')
        % axis equal        
    end % if
    
        % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/pp_error_example','.tex'], ... 
                    'dataPath',['plots/paper_plots/data/pp_error_example'], ...
                    'relativeDataPath',['figs/data/pp_error_example'], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
        
end % if

%% Alpha Gain Example

if(1)
	
    initLocalization
    
    % Plot-Calculations ___________________________________________________
    
    static_alpha = 1;
    e_g1 = 0.75;
    e_g2 = 2;
    
    e_g = [0:0.01:e_g2*1.5];
    m = -1 / (e_g2-e_g1);
    n = 1-m*e_g1;
    gain_factor = m.*e_g+n;
    gain_factor = min(max(gain_factor,0),1);    
    alpha = static_alpha .* gain_factor;
    
    % Plot ________________________________________________________________
    figure_name = ['SLERP Gain'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid off;
   
    h_plot = gobjects(0);   
    h_plot(end+1) = plot(e_g,alpha,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','alpha');    
    %axis tight       
   
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','southwest')
    xlabel('e_g [m/s^2]')
    ylabel('alpha [-]')  
    
    axes(gca)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.1, y_limits(2)+delta_y_lim*0.1];
    ylim(y_limits_new)
    vfill(e_g1)
    vfill(e_g2)
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','southwest')
    xlabel('e_g [m/s^2]')
    ylabel('alpha [-]')
   
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/slerp_gain','.tex'], ... 
                    'dataPath',['plots/paper_plots/data/slerp_gain'], ...
                    'relativeDataPath',['figs/data/slerp_gain'], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
        
end % if

%% Straight Track-Length Error Visualization

if(1)
	    
    % Plot-Calculations ___________________________________________________
    
    x_straight = 2:1:8;
    y_straight = tand(10).*x_straight;
    rng(4)
    z_data = [x_straight+1*rand(size(x_straight))-0.5; ... 
              y_straight+1*rand(size(y_straight))-0.5];
	x_straight = 0:1:10;
    y_straight = tand(10).*x_straight;
          
    
    
    % Plot ________________________________________________________________
    figure_name = ['Straight Track-Length Error Visualization'];
    close(findobj('Type','figure','Name',figure_name));
    figure('Name',figure_name); hold all; grid off;
   
    h_plot = gobjects(0);   
    h_plot(end+1) = plot(x_straight,y_straight,'k-','LineWidth',1.5,'MarkerSize',10,'DisplayName','Straight');
    h_plot(end+1) = plot(z_data(1,:),z_data(2,:),'rx','LineWidth',1.5,'MarkerSize',10,'DisplayName','z');
    %axis tight       
   
    %h_legend = legend(h_plot);
    %set(h_legend,'Location','southwest')
    xlabel('x')
    ylabel('y')
    
    axes(gca)
    x_limits = xlim; delta_x_lim = diff(x_limits);
    x_limits_new = [x_limits(1)-delta_x_lim*0.2, x_limits(2)+delta_x_lim*0.2];
    xlim(x_limits_new)
    y_limits = ylim; delta_y_lim = diff(y_limits);
    y_limits_new = [y_limits(1)-delta_y_lim*0.2, y_limits(2)+delta_y_lim*0.2];
    ylim(y_limits_new)
    
       
    
    % matlab2tikz _________________________________________________________
    if 0
        cleanfigure('targetResolution',300)
        matlab2tikz( ...
                    'filename',['plots/paper_plots/straight_length_error','.tex'], ... 
                    'dataPath',['plots/paper_plots/data/straight_length_error'], ...
                    'relativeDataPath',['figs/data/straight_length_error'], ...
                    'externalData',true, ...                
                    'floatFormat','%.8g', ... 
                    'maxChunkLength',32000, ... 
                    'width','9.5cm', ...
                    'encoding','utf8', ...
                    'standalone',false, ...
                    'figurehandle',findobj('Type','figure','Name',figure_name) ... 
                   );
    end % if
        
end % if
