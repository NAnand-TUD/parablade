function [] = plot_blade_surface(plot_options)

% Read the blade surface coordinates
file_name = 'output/coordinates/surface_coordinates.csv';
S_blade = dlmread(file_name, ',', 1, 0);
S_blade(:,[4,2]) = S_blade(:,[2,4]);

% Load the plotting settings
Nu = plot_options.Nu;
Nv = plot_options.Nv;
N_plot = plot_options.N_plot;
N_total = plot_options.N_total;
my_color = plot_options.color;

% Blade surface coordinates
[x_blade,y_blade,z_blade] = reshape_data(S_blade,Nu,Nv);

% Blade hub section coordinates
x_hub = S_blade(1:Nu,2);
y_hub = S_blade(1:Nu,3);
z_hub = S_blade(1:Nu,4);

% Blade shroud section coordinates
x_shroud = S_blade(end-Nu+1:end,2);
y_shroud = S_blade(end-Nu+1:end,3);
z_shroud = S_blade(end-Nu+1:end,4);


% Plot for axisymmetric cascades
if strcmp(plot_options.cascade_type,'axisymmetric') == 1

    % Compute the angle "alpha" spanned by the N blades
    theta = atan2(z_blade,-y_blade);
    theta_min = min(min(theta));
    theta_max = max(max(theta));
    theta_corr = (theta_min+theta_max)/2 - pi/2;
    theta_arc = 2*pi*N_plot/N_total;
    delta_theta = theta_max-theta_min;
    if N_plot == 1
        alpha = theta_corr;
    else
        alpha = linspace(theta_corr-theta_arc/2,theta_corr+theta_arc/2,N_plot);
    end
    
    % Draw a portion of the cascade by rotating the blades
    for i = 1:N_plot

        % Plot blade surface(s)
        X = x_blade;
        Y = cos(alpha(i))*y_blade - sin(alpha(i))*z_blade;
        Z = sin(alpha(i))*y_blade + cos(alpha(i))*z_blade;
        h = surf(X,Y,Z,'FaceColor',my_color,'EdgeColor','none');
        h.SpecularStrength = 0.9;
        h.SpecularExponent = 25;

        % Plot hub section(s)
        X_hub = x_hub;
        Y_hub = cos(alpha(i))*y_hub - sin(alpha(i))*z_hub;
        Z_hub = sin(alpha(i))*y_hub + cos(alpha(i))*z_hub;
        patch('XData',X_hub,'YData',Y_hub,'ZData',Z_hub,'FaceColor',my_color,'EdgeColor','none');

        % Plot shroud section(s)
        X_shroud = x_shroud;
        Y_shroud = cos(alpha(i))*y_shroud - sin(alpha(i))*z_shroud;
        Z_shroud = sin(alpha(i))*y_shroud + cos(alpha(i))*z_shroud;
        patch('XData',X_shroud,'YData',Y_shroud,'ZData',Z_shroud,'FaceColor',my_color,'EdgeColor','none');

    end
    
end


% Plot for linear cascades
if strcmp(plot_options.cascade_type,'linear') == 1

    % Compute the angle "alpha" spanned by the N blades
    y_min = min(min(y_blade));
    dy = max(max(y_blade))-min(min(y_blade));
    spacing = 1.00*dy;

    % Draw a portion of the cascade by rotating the blade
    for i = 1:N_plot
        
        % Plot blade surface(s)
        X = x_blade;
        Y = y_blade + (i-1)*spacing;
        Z = z_blade;
        surf(X,Y,Z,'FaceColor',my_color,'EdgeColor','none');
        
        % Plot hub section(s)
        X_hub = x_hub;
        Y_hub = y_hub + (i-1)*spacing;
        Z_hub = z_hub;
        patch('XData',X_hub,'YData',Y_hub,'ZData',Z_hub,'FaceColor',my_color,'EdgeColor','none');
        
        % Plot shroud section(s)
        X_shroud = x_shroud;
        Y_shroud = y_shroud + (i-1)*spacing;
        Z_shroud = z_shroud;
        patch('XData',X_shroud,'YData',Y_shroud,'ZData',Z_shroud,'FaceColor',my_color,'EdgeColor','none');

    end
    
end


end


function [X,Y,Z] = reshape_data(S,Nu,Nv)

% Preallocate space
X = zeros(Nv,Nu);
Y = zeros(Nv,Nu);
Z = zeros(Nv,Nu);

% Rename variables
x = S(:,2);
y = S(:,3);
z = S(:,4);

% Reshape the coordinates into an array
for i = 1:Nv
    X(i,:) = x(1+(i-1)*Nu:i*Nu);
    Y(i,:) = y(1+(i-1)*Nu:i*Nu);
    Z(i,:) = z(1+(i-1)*Nu:i*Nu);
end

end





