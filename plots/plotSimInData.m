% 
% Plot input data
%
%   Author: Hanno Winter
%   Date: 13-Mar-2021; Last revision: 13-Mar-2021

%% GNSS Positions

if(1)        
    
    % Pre-Calculations ____________________________________________________
        
    % Plot Routine ________________________________________________________    
    for session_i = 1
        
        % Plot-Calculations _______________________________________________
        gnss_time = datetime(seconds(gnss_data.Time),'ConvertFrom','posixtime');
        gnss_time = seconds(gnss_data.Time - gnss_data.Time(1));
        
        gnss_plot_data = gnss_data;
        numeric_columns_selector = varfun(@isnumeric,gnss_plot_data,'output','uniform');
        gnss_plot_data{~gnss_data_validity,numeric_columns_selector} = nan;
                
        min_lon = min(gnss_data.Longitude_deg);
        max_lon = max(gnss_data.Longitude_deg);
        min_lat = min(gnss_data.Latitude_deg);
        max_lat = max(gnss_data.Latitude_deg);        
               
        % Plot ____________________________________________________________        
        figure_name = ['GNSS Lat-Lon-Position (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        h_fig = figure('Name',figure_name);
        
        clear h_plot    
        h_plot = gobjects(0);
        
        if verLessThan('matlab', '9.5') || 1
            
            hold on; grid on;
            h_plot(end+1) = plot(gnss_plot_data.Longitude_deg,gnss_plot_data.Latitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
                        
            % You can try
            % https://de.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot_google_map
            % to create satellite plots with Google Maps API
            if exist('plot_google_map')
                plot_google_map('MapType', 'satellite', 'ShowLabels', 1);
            end % if    
            
            %legend(h_plot);
            xlabel('longitude [deg]')
            ylabel('latitude [deg]')
            
            dcm_obj = datacursormode(h_fig);
            set(dcm_obj,'UpdateFcn',{@myCustomDataTipFcn,{gnss_time}})
            
        else
            
            h_plot(end+1) = geoplot(gnss_plot_data.Latitude_deg,gnss_plot_data.Longitude_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS'); hold on
                        
            geobasemap('topographic');
            %legend(h_plot);
            
            dcm_obj = datacursormode(h_fig);
            set(dcm_obj,'UpdateFcn',{@myCustomDataTipFcn,{gnss_time}})
        
        end % if
    
    end % for session_i
end % if

%% GNSS Speeds

if(1)
    
    % Pre-Calculations ____________________________________________________
              
    % Plot Routine ________________________________________________________
    for session_i = 1
        
        % Plot-Calculations _______________________________________________
        gnss_time = datetime(seconds(gnss_data.Time),'ConvertFrom','posixtime');
        gnss_time = gnss_data.Time - gnss_data.Time(1);
        
        gnss_plot_data = gnss_data;
        numeric_columns_selector = varfun(@isnumeric,gnss_plot_data,'output','uniform');
        gnss_plot_data{~gnss_data_validity,numeric_columns_selector} = nan;
                
        % Plot ____________________________________________________________
        figure_name = ['GNSS NED-Velocities (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;    
        h_plot(end+1) = plot(gnss_time,gnss_plot_data.VelocityNorth_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{north} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_time,gnss_plot_data.VelocityEast_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{east} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(gnss_time,gnss_plot_data.VelocityDown_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{down} [m/s]')

        linkaxes([ax1,ax2,ax3],'x');
    end % for session_i
end % if

%% Ground Speeds (GNSS)

if(1)   
    
    % Pre-Calculations ____________________________________________________
        
    % Plot Routine ________________________________________________________
    for session_i = 1
        
        % Plot-Calculations _______________________________________________
        gnss_time = datetime(seconds(gnss_data.Time),'ConvertFrom','posixtime');
        gnss_time = gnss_data.Time - gnss_data.Time(1);
                   
        gnss_plot_data = gnss_data;
        numeric_columns_selector = varfun(@isnumeric,gnss_plot_data,'output','uniform');
        gnss_plot_data{~gnss_data_validity,numeric_columns_selector} = nan;
        
        % Plot ____________________________________________________________
        figure_name = ['GNSS Speed over Ground (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(2,1,1); hold on; grid on;
        h_plot(end+1) = plot(gnss_time,gnss_plot_data.VelocityGround_ms,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('v_{ground} [m/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(2,1,2); hold on; grid on;
        h_plot(end+1) = plot(gnss_time,gnss_plot_data.Heading_deg,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','GNSS');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('heading [deg]')

        linkaxes([ax1,ax2],'x');
    end % for session_i
end % if

%% IMU Accelerations

if(1)    
    
    % Pre-Calculations ____________________________________________________
        
    % Plot Routine ________________________________________________________
    for session_i = 1
        
        % Plot-Calculations _______________________________________________        
        imu_time = datetime(seconds(imu_data.Time),'ConvertFrom','posixtime');
        imu_time = imu_data.Time - imu_data.Time(1);
        
        % Plot ____________________________________________________________
        figure_name = ['IMU Accelerations (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(imu_time,imu_data.AccX_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMU');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('a_x [m/s^2]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(imu_time,imu_data.AccY_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMU');    
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('a_y [m/s^2]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(imu_time,imu_data.AccZ_mss,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMU');   
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('a_z [m/s^2]')

        linkaxes([ax1,ax2,ax3],'x');
    end % for session_i
end % if

%% IMU Turn Rates

if(1)  
        
    % Pre-Calculations ____________________________________________________
        
    % Plot Routine ________________________________________________________        
    for session_i = 1
        
        % Plot-Calculations _______________________________________________
        imu_time = datetime(seconds(imu_data.Time),'ConvertFrom','posixtime');
        imu_time = imu_data.Time - imu_data.Time(1);
        
        % Plot ____________________________________________________________
        figure_name = ['IMU Turn Rates (Session ',sprintf('%02i',session_i),')'];
        close(findobj('Type','figure','Name',figure_name));
        figure('Name',figure_name); hold all; grid on;

        clear h_plot    
        h_plot = gobjects(0);   
        ax1 = subplot(3,1,1); hold on; grid on;
        h_plot(end+1) = plot(imu_time,imu_data.TurnRateX_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMU');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('w_x [deg/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax2 = subplot(3,1,2); hold on; grid on;
        h_plot(end+1) = plot(imu_time,imu_data.TurnRateY_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMU');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('w_y [deg/s]')

        clear h_plot
        h_plot = gobjects(0); 
        ax3 = subplot(3,1,3); hold on; grid on;
        h_plot(end+1) = plot(imu_time,imu_data.TurnRateZ_degs,'.-','LineWidth',1.5,'MarkerSize',10,'DisplayName','IMU');
        %h_legend = legend(h_plot);
        %set(h_legend,'Location','southwest')
        xlabel('UTC time')
        ylabel('w_z [deg/s]')

        linkaxes([ax1,ax2,ax3],'x');
    end % for session_i
end % if

%% Helper Functions

function txt = myCustomDataTipFcn(pointDataTip,event_obj,utc_time)

data_index = pointDataTip.Cursor.DataIndex;
pos = get(event_obj,'Position');

% time_str = datestr(utc_time{1}(data_index));
time_str = num2str(utc_time{1}(data_index));

txt = { ...
        ['X: ',num2str(pos(1))], ...
        ['Y: ',num2str(pos(2))],  ... 
        ['Time: ',time_str] ... 
      };

end 
