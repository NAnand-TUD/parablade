function [] = plot_blade_normals(scale)

% Plot the unitary vectors normal to the blade surface as straigth 
% lines whose size is controlled by the input 'scale'

% Load the surface coordinates
file_name = 'output/sensitivities/grad_stagger_0.csv';
S_blade = dlmread(file_name, ',', 1, 0);
C = S_blade(:,2:4);
N = S_blade(:,8:10);

% Arrange the order for SU2 convention
C(:,[1,3]) = C(:,[3,1]);
N(:,[1,3]) = N(:,[3,1]);
N = N*scale;

% Plot the normal vectors
for i = 1:length(S_blade(:,1))

    X = [C(i,1) C(i,1) + N(i,1)];
    Y = [C(i,2) C(i,2) + N(i,2)];
    Z = [C(i,3) C(i,3) + N(i,3)];

    plot3(X,Y,Z, ...
        'Marker', 'none', ...
        'Color', 'k', ...
        'LineStyle', '-', ...
        'MarkerSize', 1.00, ...
        'MarkerFaceColor', 'b', ...
        'LineWidth', 0.50);
        

    
end
