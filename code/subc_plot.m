clear; clc; close all;
addpath results\

% Specify the path to your CSV file
csv_file_path = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\atlas_plot\data\plot_subcortical_data.csv';

% Read the CSV file into MATLAB
data_table = readtable(csv_file_path);

% Create a variable names array for iteration
variables = {'t', 'CID', 'HCs', 'Sig','RA_T', 'RA_Pre', 'RA_Pos', 'RA_Sig','SA_T', 'SA_Pre', 'SA_Pos', 'SA_Sig'};

% Define fixed color ranges for each variable outside the loop
color_ranges = {
    [-4.0 4.0],  % Color range for 't'
    [38.0 59.0],  % Color range for 'CID'
    [38.0 59.0],  % Color range for 'HCs'
    [-4.0 4.0],  % Color range for 't'
     % Color range for 't'
    [-2.8 2.8],  % Color range for 't'
    [38.0 59.0],  % Color range for 'CID'
    [38.0 59.0],  % Color range for 'HCs'
    [-2.8 2.8] % Color range for 'Sig'
    % Color range for 'HCs'
    [-2.8 2.8],  % Color range for 't'
    [38.0 59.0],  % Color range for 'CID'
    [38.0 59.0],  % Color range for 'HCs'
    [-2.8 2.8]
};

% Define fixed colormaps for each variable outside the loop
colormaps = {
    'RdBu_r',    % Colormap for 't'
    'autumn_r',  % Colormap for 'CID'
    'autumn_r',  % Colormap for 'HCs'
    'RdBu_r',     % Colormap for 'Sig'
    % Colormap for 'Sig'
    'RdBu_r',    % Colormap for 't'
    'autumn_r',  % Colormap for 'CID'
    'autumn_r',  % Colormap for 'HCs'
    'RdBu_r',    % Colormap for 'Sig'
        % Colormap for 'Sig'
    'RdBu_r',    % Colormap for 't'
    'autumn_r',  % Colormap for 'CID'
    'autumn_r',  % Colormap for 'HCs'
    'RdBu_r'     % Colormap for 'Sig'
};

% Define a base filename for saving the plots
base_filename = 'results\subcortical_';

% Loop through each variable
for i = 1:length(variables)
    % Select the data for the current variable
    current_data = data_table.(variables{i});
    
    % Get the color range and colormap for the current variable
    color_range = color_ranges{i};
    cmap = colormaps{i};
    
    % Plot subcortical data for the current variable
    f1 = figure('Renderer', 'painters', 'Position', [10 10 1000 667]);
    plot_subcortical(current_data, 'color_range', color_range, 'cmap', cmap, 'ventricles', 'False');
    colorbar('off'); % Turn off the colorbar display
    drawnow;
    ax = gca;
    ax.FontSize = 16; % Change the font size of the axis here

    % Construct the filename for the current plot
    filename1 = [base_filename, variables{i}, '.pdf'];
    filename2 = [base_filename, variables{i}, '.png'];
    % Save the current plot as a PDF
    exportgraphics(f1, filename1, 'ContentType', 'vector', 'Resolution', 300);
    % Save the current plot as a PNG
    exportgraphics(f1, filename2, 'Resolution', 600);  % PNG format with specified resolution
       
    % Close the figure window to free up system resources
    close(f1);
end
