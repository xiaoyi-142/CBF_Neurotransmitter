clear;clc;close all;
% 从 Excel 文件中读取数据
filePath = 'sig_region&scales_data.xlsx';
data = readtable(filePath,'Sheet', 1);

% 获取显著区域的名称
significantRegions = data.Properties.VariableNames(2:end-2); % 最后两列假设是 PSQI 和 ISI


% 分别提取 PSQI 和 ISI 数据
psqiData = data.PSQI;
isiData = data.ISI;

% 初始化存储原始 p 值的数组
p_values_psqi = zeros(length(significantRegions), 1);
p_values_isi = zeros(length(significantRegions), 1);
r_values_psqi = zeros(length(significantRegions), 1);
r_values_isi = zeros(length(significantRegions), 1);

% 计算每个显著区域的相关性和 p 值
for i = 1:length(significantRegions)
    cbfValues = data{:, significantRegions{i}};
    [r_values_psqi(i), p_values_psqi(i)] = corr(cbfValues, psqiData, 'Rows', 'complete');
    [r_values_isi(i), p_values_isi(i)] = corr(cbfValues, isiData, 'Rows', 'complete');
end

% FDR 矫正
p_fdr_psqi = mafdr(p_values_psqi, 'BHFDR', true);
p_fdr_isi = mafdr(p_values_isi, 'BHFDR', true);

% 对 PSQI 进行分析
for i = 1:length(significantRegions)
    cbfValues = data{:, significantRegions{i}};
    significantRegions_name{i} = strrep(significantRegions{i}, '_', ' ');
    % 检查 PSQI 的 FDR 矫正后的 p 值是否显著
    if p_fdr_psqi(i) < 0.05
        figure;
        scatter(psqiData, cbfValues, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7]);

        % 回归线
        hold on; % 保持当前图形，用于添加更多图层
        coeffs = polyfit(psqiData, cbfValues, 1);
        fittedX = linspace(min(psqiData), max(psqiData), 200);
        fittedY = polyval(coeffs, fittedX);
        plot(fittedX, fittedY, 'Color', [1, 0.3, 0], 'LineWidth', 1.5);
        hold off;
        % 标注相关性和 P 值
        text_pos_x = min(psqiData) + 0.25 * range(psqiData);
        text_pos_y = min(cbfValues) + 0.10 * range(cbfValues);
        corrText = sprintf('r = %.2f, $p_{\\mathrm{FDR}} = %.3f$', r_values_psqi(i),  p_fdr_psqi(i));
 
        text(text_pos_x, text_pos_y, corrText, 'FontSize', 18, 'FontName', 'Arial', 'HorizontalAlignment', 'center','Interpreter', 'latex');

        %         text(text_pos_x, text_pos_y, sprintf('r = %.2f, pfdr = %.3f', r_values_psqi(i), p_fdr_psqi(i)), 'FontSize', 18);
        set(gca, 'FontSize', 18);

        % 标签和标题
        xlabel('PSQI scores', 'FontSize', 18);
        ylabel([significantRegions_name{i} ' CBF values'], 'FontSize', 18);

        % 保存图像
        print('-dpng', fullfile('scatter_plots', ['PSQI_', strrep(significantRegions_name{i},' ', '_'), '.png']), '-r600');
        close;
    end
end

% 对 ISI 进行分析
for i = 1:length(significantRegions)
    cbfValues = data{:, significantRegions{i}};
    significantRegions_name{i} = strrep(significantRegions{i}, '_', ' ');
    % 检查 ISI 的 FDR 矫正后的 p 值是否显著
    if p_fdr_isi(i) < 0.05
        figure;
        scatter(isiData, cbfValues, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7]);

        % 回归线
        hold on;
        coeffs = polyfit(isiData, cbfValues, 1);
        fittedX = linspace(min(isiData), max(isiData), 200);
        fittedY = polyval(coeffs, fittedX);
        plot(fittedX, fittedY, 'Color', [1, 0.3, 0], 'LineWidth', 1.5);
        hold off;
        % 标注相关性和 P 值
        text_pos_x = min(isiData) + 0.4 * range(isiData);
        text_pos_y = min(cbfValues) + 0.02 * range(cbfValues);
        corrText = sprintf('r = %.2f, $p_{\\mathrm{FDR}} = %.3f$', r_values_isi(i),  p_fdr_isi(i));
        text(text_pos_x, text_pos_y, corrText, 'Interpreter', 'latex');

%         text(text_pos_x, text_pos_y, sprintf('r = %.2f, $p_{\\mathrm{fdr}} = %.3f', r_values_isi(i), p_fdr_isi(i)), 'FontSize', 18, 'Interpreter', 'latex');
        set(gca, 'FontSize', 18);
        % 标签和标题
        xlabel('ISI scores', 'FontSize', 18);
        ylabel([significantRegions_name{i} ' CBF values'], 'FontSize', 18);

        % 保存图像
        print('-dpng', fullfile('scatter_plots', ['ISI_', strrep(significantRegions_name{i},' ', '_'), '.png']), '-r300');
        close;
    end
end



%%对量表的相关性进行绘制
% 初始化相关性系数和 p 值数组
numRegions = length(significantRegions);
corrCoeffs = zeros(numRegions, 1);
pValues = ones(numRegions, 1);

% 计算 PSQI 与每个区域的相关性和 p 值
for i = 1:numRegions
    [r, p] = corr(psqiData, data{:, significantRegions{i}}, 'Rows', 'complete');
    corrCoeffs(i) = r;
    pValues(i) = p;
end

% 进行 FDR 矫正
pValuesFDR = mafdr(pValues, 'BHFDR', true);

% 筛选显著的相关性
significantIdx = pValuesFDR < 0.01;
significantIdx = {true, true, true, false, true, false, false, false, false, false, false, false, false, false, false, false};
significantIdx = cell2mat(significantIdx);

significantCorrCoeffs = corrCoeffs(significantIdx);
significantRegionsFiltered = significantRegions(significantIdx);
% significantRegionsFiltered 中的原始名称
% significantRegionsFiltered = {'PUTAMEN_LEFT', 'DIENCEPHALON_VENTRAL_LEFT', 'PALLIDUM_RIGHT', 'DIENCEPHALON_VENTRAL_RIGHT', 'RH_Vis_4'};

% 为这些区域赋予自定义的名称
customNames = {'L Put', 'L Pal', 'L DEV', 'R DEV'};
customNames = {'L DEV', 'R DEV', 'L Pal', 'L Put'};
% 确保 significantCorrCoeffs 至少有一个元素
if isempty(significantCorrCoeffs)
    error('没有显著的相关性系数。');
end

% 确保 customNames 数组的长度与 significantCorrCoeffs 一致
if length(customNames) ~= length(significantCorrCoeffs)
    error('customNames 的长度必须与 significantCorrCoeffs 的长度相同。');
end

% 将 significantCorrCoeffs 转换为行向量
significantCorrCoeffsRow = significantCorrCoeffs';

% 绘制雷达图
figure;
pax = polaraxes;
hold on;

% 定义颜色和样式
lineColor = [0, 0.4470, 0.7410]; % 蓝色
markerColor = [0.8500, 0.3250, 0.0980]; % 橙色
fillColor = [0.9290, 0.6940, 0.1250]; % 黄色

% 极坐标转换为笛卡尔坐标并绘制雷达图
theta = linspace(0, 2*pi, numel(significantCorrCoeffsRow)+1);
rho = [significantCorrCoeffsRow, significantCorrCoeffsRow(1)];
[X, Y] = pol2cart(theta, rho);
polarplot(pax, theta, rho, 'o-', 'LineWidth', 2, 'Color', lineColor, 'MarkerFaceColor', markerColor);

% 绘制连接的线段来形成填充区域
for i = 1:length(X)-1
    polarplot(pax, [theta(i), theta(i+1)], [rho(i), rho(i+1)], 'LineWidth', 2, 'Color', fillColor);
end

% 填充最后一段线
polarplot(pax, [theta(end), theta(1)], [rho(end), rho(1)], 'LineWidth', 2, 'Color', fillColor);

% 添加显著性标记
for i = 1:length(X)-1
    if p_fdr_psqi(i) < 0.1 % 显著性检验
        text(theta(i), rho(i)-0.03, '*', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 20);
    end
end

% 设置雷达图属性
pax.ThetaTick = linspace(0, 360, numel(customNames) + 1);
pax.ThetaTickLabel = [customNames, customNames{1}];
pax.RTick = 0.1:0.1:0.36; % 根据需要调整
pax.RLim = [0, 0.36];   % 根据需要调整

% 图表标题
title('PSQI negative correlation values', 'FontSize', 22);

% 可视化优化
set(pax, 'FontSize', 18);
set(pax, 'ThetaZeroLocation', 'top');

% 调整图表比例和布局
set(pax, 'GridLineStyle', '--');
set(pax, 'GridAlpha', 0.5); % 调整网格线透明度

% 输出雷达图到指定位置
outputFilePath = 'scatter_plots/render_plot.png'; % 替换为您的文件路径
exportgraphics(pax, outputFilePath, 'Resolution', 600);
