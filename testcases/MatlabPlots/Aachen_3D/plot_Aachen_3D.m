%% ParaBlade/Matlab plots

% Clean workspace
clear all
close all
clc

% Add plotting functions
addpath(genpath('../../../src/src_matlab_plot'))
plot_settings

% Prepare directory to save figures
if ~exist('MATLAB_figures', 'dir')
   mkdir MATLAB_figures
end

N_blades = 50;
Nu = 500;       % hardcoded, see Blade3D initialize_uv_values()
Nv = 25;        % hardcoded, see Blade3D initialize_uv_values()
cascade_type = 'axisymmetric';


%% Plot the surface of the blade
% Prepare the figure
figure1 = figure(); ax_fig1 = gca;
hold on; axis image; box on;
% view(300, 25); camlight('headlight'); lighting flat; view(-85,40)
view(-20, -10); camlight('headlight'); lighting flat; view(-75,30)
xlabel('$x$ axis');
ylabel('$y$ axis');
zlabel('$z$ axis');
ax_fig1.XAxis.TickLabelFormat = '%.2f';
ax_fig1.YAxis.TickLabelFormat = '%.2f';
ax_fig1.ZAxis.TickLabelFormat = '%.2f';
axis off

% Prepare plotting options
plot_options.cascade_type = cascade_type;
plot_options.Nu = Nu;
plot_options.Nv = Nv;
plot_options.N_total = N_blades;
plot_options.N_plot = 5;
plot_options.color = 0.6*[1 1 1];

% Make the plot
plot_blade_surface(plot_options)
plot_hub_surface(plot_options)
% plot_shroud_surface(plot_options)

% % Save the figure
% export_fig(figure1,'./MATLAB_figures/Aachen_3D_b','-png','-r1000')



%% Plot the surface sensitivity
% % Prepare the figure
% figure2 = figure(); ax_fig2 = gca;
% hold on; axis image; box on; view(290,30);
% % camlight('headlight'); lighting flat;
% xlabel('$x$ axis');
% ylabel('$y$ axis');
% zlabel('$z$ axis');
% ax_fig2.XAxis.TickLabelFormat = '%.2f';
% ax_fig2.YAxis.TickLabelFormat = '%.2f';
% ax_fig2.ZAxis.TickLabelFormat = '%.2f';
% my_cmap = parula(100);
% colormap(my_cmap)
% shading interp
% % shading flat
% axis off
% 
% % Prepare plotting options
% variable_name = 'dist_in_0';
% plot_options.N_sections = N_sections;
% plot_options.N_total = N_blades;
% plot_options.N_plot = 1;
% plot_options.color = my_cmap(floor(end/2),:);
% plot_options.colorbar = 'no';
% 
% % Make the plot
% plot_blade_sensitivity(variable_name, plot_options)
% % plot_hub_surface(plot_options)
% % plot_shroud_surface(plot_options)
% 
% % Save the figure
% % export_fig(figure2,'./matlab_figures/NASA_R67_sensitivity','-png','-r1000')


