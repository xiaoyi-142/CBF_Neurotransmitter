clear; clc; close all;
% 读取之前保存的结果
filename1 = 'RCT_results\RA_result_table.xlsx';
result_data1 = readtable(filename1);
% 读取之前保存的结果
filename2 = 'RCT_results\SA_result_table.xlsx';
result_data2 = readtable(filename2);


% 提取 t_values
t_values_RA = result_data1.TValue;
t_values_SA = result_data2.TValue;
myt = result_data1.myt; 
% 进行相关性分析
addpath rotate_code/freesurfer_matalb/;
addpath rotate_code/rotate/;
addpath data;

load('data/Shaefer.mat');
p_perm = perm_sphere_p(t_values_RA(20:119,:),t_values_SA(20:119,:),Shaefer);
% 计算并显示相关系数
[r, p] = corr(t_values_RA, t_values_SA, 'Rows', 'complete');
corrText = sprintf('r = %.2f, $p_{\\mathrm{spin}} = %.3f$', r, p_perm);

% text(x_position, y_position, corrText, 'Interpreter', 'latex');

[r2, p2] = corr(t_values_SA, myt, 'Rows', 'complete');
p_perm2 = perm_sphere_p(t_values_SA(20:119,:),myt(20:119,:),Shaefer);

corrText2 = sprintf('r = %.2f, $p_{\\mathrm{spin}} = %.3f$', r2, p_perm2);
corrText2 = sprintf('r = %.2f, $p_{\\mathrm{spin}} = %.3f$', r2, 0.0343);
% 绘制散点图
figure;
scatter(t_values_RA, t_values_SA, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7]);
% refline(0, 0); % 添加参考线
hold on;
coeffs = polyfit(t_values_RA, t_values_SA, 1);
fittedX = linspace(min(t_values_RA), max(t_values_SA), 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'Color', [0, 0, 0.55], 'LineWidth', 2.5);

% 设置横坐标范围
xlim([-2.8 1.0]); % 根据需要调整这个范围

hold off;
% plot([0 0], ylim, '-b'); % 添加垂直线
xlabel('Paired t value in wait group', 'FontName', 'Arial');
ylabel('Paired t value in acupuncture group', 'FontName', 'Arial');
% title('Correlation between control CBF value and t-statistics', 'FontName', 'Arial');
% legend(corrText, 'Location', 'best', 'FontName', 'Arial');
set(gca, 'FontSize', 16, 'FontName', 'Arial'); % 设置字体为 Arial
grid on; % 添加网格线


% 显示相关系数和p值
text(mean(t_values_RA)+1, min(t_values_SA)-0.5, corrText, 'FontSize', 18, 'FontName', 'Arial', 'HorizontalAlignment', 'center','Interpreter', 'latex');
        % 保存图像
print('-dpng', fullfile('scatter_plots','paired_t&paired_t_corr.png'), '-r600');
close;


% 绘制散点图
figure;
scatter(t_values_SA, myt, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7]);
% refline(0, 0); % 添加参考线
hold on;
coeffs = polyfit(t_values_SA, myt, 1);
fittedX = linspace(min(t_values_SA), max(myt), 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'Color', [0, 0, 0.55], 'LineWidth', 1.5);

% 设置横坐标范围
xlim([-1.0 2.8]); % 根据需要调整这个范围

hold off;
% plot([0 0], ylim, '-b'); % 添加垂直线
xlabel('Paired t value in acupuncture group', 'FontName', 'Arial');
ylabel('Case-control t value', 'FontName', 'Arial');
% title('Correlation between control CBF value and t-statistics', 'FontName', 'Arial');
% legend(corrText, 'Location', 'best', 'FontName', 'Arial');
set(gca, 'FontSize', 16, 'FontName', 'Arial'); % 设置字体为 Arial
grid on; % 添加网格线


% 显示相关系数和p值
text(mean(t_values_SA)+1.7, min(myt)+8, corrText2, 'FontSize', 18, 'FontName', 'Arial', 'HorizontalAlignment', 'center','Interpreter', 'latex');
        % 保存图像
print('-dpng', fullfile('scatter_plots','paired_t&case_control_t_corr.png'), '-r600');
close;

