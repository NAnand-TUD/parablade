function [] = plot_blade_coordinates()

% Plots the coordinates of the blade surface as points
file_name = 'output/coordinates/surface_coordinates.csv';
S = dlmread(file_name, ',', 1, 0);

% Plot the blade coordinates as an scattered plot
plot3(S(:,2),S(:,3),S(:,4), ...
    'Marker', 'o', ...
    'Color', 'k', ...
    'LineStyle', 'none', ...
    'MarkerSize', 1.00, ...
    'MarkerFaceColor', 'k');


end

