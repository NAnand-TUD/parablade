function [] = plot_hub_surface(plot_options)

% Read the hub surface coordinates
file_name = 'output/coordinates/hub_coordinates.csv';
S_hub = dlmread(file_name, ',', 1, 0);
S_hub(:,[4,2]) = S_hub(:,[2,4]);
[N_hub,~] = size(S_hub);

% Load the plotting settings
N_plot = plot_options.N_plot;
N_total = plot_options.N_total;
my_color = plot_options.color;

% Hub section coordinates
x_hub = S_hub(1:N_hub,2);
y_hub = S_hub(1:N_hub,3);
z_hub = S_hub(1:N_hub,4);

% Plot for axisymmetric cascades
if strcmp(plot_options.cascade_type,'axisymmetric') == 1

    % Plot the hub surface
    N_theta = 250;
    theta = get_theta_angle(N_plot, N_total, N_theta);
    X_hub = zeros(N_theta,N_hub);
    Y_hub = zeros(N_theta,N_hub);
    Z_hub = zeros(N_theta,N_hub);
    
    for i = 1:N_theta
        X_hub(i,:) = x_hub;
        Y_hub(i,:) = cos(theta(i))*y_hub - sin(theta(i))*z_hub;
        Z_hub(i,:) = sin(theta(i))*y_hub + cos(theta(i))*z_hub;
    end
    
    surf(X_hub,Y_hub,Z_hub,'FaceColor',my_color,'EdgeColor','none');

end


% Plot for linear cascades
if strcmp(plot_options.cascade_type,'linear') == 1
    
    % Plot the hub surface
    Ny = 2;
    y_hub = get_y_length(N_plot, Ny);
    
    X_hub = zeros(Ny,N_hub);
    Y_hub = zeros(Ny,N_hub);
    Z_hub = zeros(Ny,N_hub);
    
    for i = 1:Ny
        X_hub(i,:) = x_hub;
        Y_hub(i,:) = y_hub(i);
        Z_hub(i,:) = z_hub;
    end
    
    surf(X_hub,Y_hub,Z_hub,'FaceColor',my_color,'EdgeColor','none');

end

end


function theta = get_theta_angle(N_plot, N_total, N_theta)

% Read the blade surface coordinates
file_name = 'output/coordinates/surface_coordinates.csv';
S_blade = dlmread(file_name, ',', 1, 0);
z_blade = S_blade(:,2);
y_blade = S_blade(:,3);

% Compute the angle "theta" to draw the hub/shroud surfaces
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


















