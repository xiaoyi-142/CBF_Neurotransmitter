clear; clc; close all;

% Define the sets
case_control = {'M1', 'MOR', 'GABAa', 'VAChT', 'H3', '5HTT', '5HT1a', 'D1', 'CB1'};
paired_t = {'GABAa', 'α4β2', 'D1', '5HT4', 'mGluR5', '5HT1b'};
% case_control = {'L Put', 'L Pal', 'L DEV', 'R Pal', 'R DEV', 'L Vis 4', 'L Vis 5', 'L Vis 8', 'L Vis 9', 'L Default PFC', 'R Vis 4', 'R Vis 5', 'R Vis 6', 'R DorsAttn Post 1', 'R DorsAttn Post 4', 'R Cont Par'};
% case_control = { };
% paired_t = { };
% Find common elements
common = intersect(case_control, paired_t);

% Find unique elements
unique_case_control = setdiff(case_control, paired_t);
unique_paired_t = setdiff(paired_t, case_control);

% Create a figure
figure;
hold on;

% Plot circles for Venn diagram
theta = linspace(0, 2*pi, 100);
x1 = 0.5 + 0.5 * cos(theta);
y1 = 0.5 + 0.5 * sin(theta);
x2 = 1.0 + 0.5 * cos(theta);
y2 = 0.5 + 0.5 * sin(theta);

fill(x1, y1, 'b', 'FaceAlpha', 0.5); % case-control t circle
fill(x2, y2, 'r', 'FaceAlpha', 0.5); % paired t circle
% fill(x1, y1, [0.3010, 0.7450, 0.9330], 'FaceAlpha', 0.5); % case-control t circle (Light Blue)
% fill(x2, y2, [1, 0.4, 0.4], 'FaceAlpha', 0.5); % paired t circle (Light Red)


% Define font parameters
font_name = 'Arial';
font_size = 14;

% Add labels for unique and common elements
text(0.3, 0.5, sprintf('%s\n', unique_case_control{:}), 'Color', 'w', 'FontName', font_name, 'FontSize', font_size);
text(1.2, 0.5, sprintf('%s\n', unique_paired_t{:}), 'Color', 'w', 'FontName', font_name, 'FontSize', font_size);
text(0.75, 0.5, sprintf('%s\n', common{:}), 'Color', 'w', 'FontName', font_name, 'FontSize', font_size);

% Add subtitles for each circle outside the circles
text(0.2, -1.2, 'Case-Control t', 'Color', 'k', 'FontWeight', 'bold', 'FontName', font_name, 'FontSize', font_size);
text(0.6, -1.2, 'Paired t', 'Color', 'k', 'FontWeight', 'bold', 'FontName', font_name, 'FontSize', font_size);

% Set plot limits and turn off axes
axis equal;
axis off;

% Add title
title('Neurotransmitters contributors in both case-Control t and paired t', 'FontName', font_name, 'FontSize', font_size);

hold off;

% Save the figure
if ~exist('scatter_plots', 'dir')
    mkdir('scatter_plots')
end
% print('scatter_plots/venn_diagram_neurotransmitter', '-dpng', '-r600');
print('scatter_plots/venn_diagram_region', '-dpng', '-r600');