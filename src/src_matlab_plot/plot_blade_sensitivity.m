function [] = plot_blade_sensitivity(variable_name, plot_options)

% Read the surface coordinates and sensitivity
file_name = ['output/sensitivities/grad_',variable_name,'.csv'];
S_blade = dlmread(file_name, ',', 1, 0);

% Arrange the order for SU2 convention
S_blade(:,[4,2]) = S_blade(:,[2,4]);
S_blade(:,[7,5]) = S_blade(:,[5,7]);
S_blade(:,[10,8]) = S_blade(:,[8,10]);

% Load the plotting settings
Nu = plot_options.Nu;
Nv = plot_options.Nv;
N_plot = plot_options.N_plot;
N_total = plot_options.N_total;
my_color = plot_options.color;

% Blade surface coordinates and surface sensitivity
sens = dot(S_blade(:,5:7), S_blade(:,8:10), 2);
[x_blade,y_blade,z_blade,sens] = reshape_data([S_blade(:,2:4),sens],Nu,Nv);

% Hub section coordinates
x_hub = S_blade(1:Nu,2);
y_hub = S_blade(1:Nu,3);
z_hub = S_blade(1:Nu,4);

% Shroud section coordinates
x_shroud = S_blade(end-Nu+1:end,2);
y_shroud = S_blade(end-Nu+1:end,3);
z_shroud = S_blade(end-Nu+1:end,4);

% Adjust colormap options
a = min(min(sens));
b = max(max(sens));
colormap_bound = max(abs(a),abs(b));
caxis([-colormap_bound,colormap_bound]);

if strcmp(plot_options.colorbar,'yes') == 1
    
    % Create a color bar
    c = colorbar;
    c.Label.Interpreter = 'Latex';
    c.TickLabelInterpreter = 'Latex';
    c.FontSize = 10;

    % Change the number of decimals of the colormap legend
    num = c.Ticks;
    cell_str = cell(length(num));
    for i = 1:length(num)
        cell_str{i} = num2str(num(i),'%0.2f');
    end
    c.TickLabels = cell_str;
    
end


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

    % Draw a portion of the cascade by rotating the blade
    for i = 1:N_plot

        % Plot blade surface(s)
        X = x_blade;
        Y = cos(alpha(i))*y_blade - sin(alpha(i))*z_blade;
        Z = sin(alpha(i))*y_blade + cos(alpha(i))*z_blade;
        surf(X,Y,Z,sens,'EdgeColor','none');

        % Plot hub section(s)
        X_hub = x_hub;
        Y_hub = cos(alpha(i))*y_hub - sin(alpha(i))*z_hub;
        Z_hub = sin(alpha(i))*y_hub + cos(alpha(i))*z_hub;
        patch('XData',X_hub,'YData',Y_hub,'ZData',Z_hub,'FaceColor',my_color,'EdgeColor','k');

        % Plot shroud section(s)
        X_shroud = x_shroud;
        Y_shroud = cos(alpha(i))*y_shroud - sin(alpha(i))*z_shroud;
        Z_shroud = sin(alpha(i))*y_shroud + cos(alpha(i))*z_shroud;
        patch('XData',X_shroud,'YData',Y_shroud,'ZData',Z_shroud,'FaceColor',my_color,'EdgeColor','k');

    end
end


% Plot for linear cascades
if strcmp(plot_options.cascade_type,'linear') == 1

    % Compute the angle "alpha" spanned by the N blades
    dx = max(max(x_blade))-min(min(x_blade));
    spacing = 0.75*dx;

    % Draw a portion of the cascade by rotating the blade
    for i = 1:N_plot
        % Blade surface
        X = x_blade;
        Y = y_blade + (i-1)*spacing;
        Z = z_blade;
        surf(X,Y,Z,sens,'EdgeColor','none');
        
        X_hub = x_hub;
        Y_hub = y_hub + (i-1)*spacing;
        Z_hub = z_hub;
        patch('XData',X_hub,'YData',Y_hub,'ZData',Z_hub,'FaceColor',my_color,'EdgeColor','k');
        
        X_shroud = x_shroud;
        Y_shroud = y_shroud + (i-1)*spacing;
        Z_shroud = z_shroud;
        patch('XData',X_shroud,'YData',Y_shroud,'ZData',Z_shroud,'FaceColor',my_color,'EdgeColor','k');

    end
    
end


end


function [X,Y,Z,G] = reshape_data(S,Nu,Nv)

% Preallocate space
X = zeros(Nv,Nu);
Y = zeros(Nv,Nu);
Z = zeros(Nv,Nu);
G = zeros(Nv,Nu);

% Rename variables
x = S(:,1);
y = S(:,2);
z = S(:,3);
g = S(:,4);

% Reshape the coordinates into an array
for i = 1:Nv
    X(i,:) = x(1+(i-1)*Nu:i*Nu);
    Y(i,:) = y(1+(i-1)*Nu:i*Nu);
    Z(i,:) = z(1+(i-1)*Nu:i*Nu);
    G(i,:) = g(1+(i-1)*Nu:i*Nu);
end

end

