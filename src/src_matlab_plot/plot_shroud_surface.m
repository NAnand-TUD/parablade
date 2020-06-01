function [] = plot_shroud_surface(plot_options)

% Read the shroud surface coordinates
file_name = 'output/coordinates/shroud_coordinates.csv';
S_shroud = dlmread(file_name, ',', 1, 0);
S_shroud(:,[4,2]) = S_shroud(:,[2,4]);
[N_shroud,~] = size(S_shroud);

% Load the plotting settings
N_plot = plot_options.N_plot;
N_total = plot_options.N_total;
my_color = plot_options.color;

% shroud section coordinates
x_shroud = S_shroud(1:N_shroud,2);
y_shroud = S_shroud(1:N_shroud,3);
z_shroud = S_shroud(1:N_shroud,4);

% Plot for axisymmetric cascades
if strcmp(plot_options.cascade_type,'axisymmetric') == 1
        
    % Plot the shroud surface
    N_theta = 250;
    theta = get_theta_angle(N_plot, N_total, N_theta);
    X_shroud = zeros(N_theta,N_shroud);
    Y_shroud = zeros(N_theta,N_shroud);
    Z_shroud = zeros(N_theta,N_shroud);
    
    for i = 1:N_theta
        X_shroud(i,:) = x_shroud;
        Y_shroud(i,:) = cos(theta(i))*y_shroud - sin(theta(i))*z_shroud;
        Z_shroud(i,:) = sin(theta(i))*y_shroud + cos(theta(i))*z_shroud;
    end
    
    surf(X_shroud,Y_shroud,Z_shroud,'FaceColor',my_color,'EdgeColor','none');

end


% Plot for linear cascades
if strcmp(plot_options.cascade_type,'linear') == 1
    
    % Plot the shroud surface
    Ny = 2;
    y_shroud = get_y_length(N_plot, Ny);
    
    X_shroud = zeros(Ny,N_shroud);
    Y_shroud = zeros(Ny,N_shroud);
    Z_shroud = zeros(Ny,N_shroud);
    
    for i = 1:Ny
        X_shroud(i,:) = x_shroud;
        Y_shroud(i,:) = y_shroud(i);
        Z_shroud(i,:) = z_shroud;
    end
    
    surf(X_shroud,Y_shroud,Z_shroud,'FaceColor',my_color,'EdgeColor','none');

end

end


function theta = get_theta_angle(N_plot, N_total, N_theta)

% Read the blade surface coordinates
file_name = 'output/coordinates/surface_coordinates.csv';
S_blade = dlmread(file_name, ',', 1, 0);
z_blade = S_blade(:,2);
y_blade = S_blade(:,3);

% Compute the angle "theta" to draw the shroud/shroud surfaces
theta = atan2(z_blade,-y_blade);
theta_min = min(min(theta));
theta_max = max(max(theta));
theta_arc = 2*pi*N_plot/N_total;
delta_theta = theta_max-theta_min;
a = 0.75*delta_theta;
theta = linspace(-theta_arc/2-a,theta_arc/2+a,N_theta);

end


function y = get_y_length(N_plot, Ny)

% Read the blade surface coordinates
file_name = 'output/coordinates/surface_coordinates.csv';
S_blade = dlmread(file_name, ',', 1, 0);
y_blade = S_blade(:,3);

% Compute the length "y" to draw the hub/shroud surfaces
y_min = min(min(y_blade));
dy = max(max(y_blade))-min(min(y_blade));
spacing = 1.00*dy;
y = linspace(y_min-dy/2, y_min+(N_plot-1)*spacing+5/4*dy,Ny);

end




