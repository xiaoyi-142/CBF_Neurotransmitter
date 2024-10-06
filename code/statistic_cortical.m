clear;clc;close all;
% 指定Excel文件路径
% filePath = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\HC&CID_scheaf100.xlsx';
filePath = "C:\Users\FEARLESS\Desktop\CID_ASL\source_data\Glasser_data\HC&CID_Glasser360.xlsx";
% 读取Excel数据
data1 = readtable(filePath, 'Sheet', 1); % 读取第一张工作表
data2 = readtable(filePath, 'Sheet', 2); % 读取第二张工作表
data3 = readtable(filePath, 'Sheet', 4); % 读取第三张工作表
data4 = readtable(filePath, 'Sheet', 5); % 读取第三张工作表
% 合并两个数据表，然后将结果与第三个表合并
mergedDataTemp = join(data1, data2, 'Keys', 'SubjectID');
mergedData = join(mergedDataTemp, data3, 'Keys', 'SubjectID');

% 过滤掉包含缺失值的记录
mergedData = rmmissing(mergedData);



% 初始化年龄、性别、组别和CBF数据
age = mergedData.Age;
sex = mergedData.Gender;
group = mergedData.Group;
cbfData_raw = mergedData{:, 2:380}; % 假设CBF数据是从第二列到第120列
cbfData = zscore(cbfData_raw);
% 初始化结果变量
nregs = size(cbfData, 2); % 脑区数量
mytstat = zeros(nregs, 1);
mypval = zeros(nregs, 1);

% 线性回归分析
for region = 1:nregs    
    tbl = table(age, sex, group, cbfData(:, region), 'VariableNames', {'age', 'sex', 'group', 'CBF'});
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    % 假设 tbl.group 已经是一个分类变量
    tbl.group = categorical(tbl.group, {'HC', 'CID'}); % 确保'HC'是第一个类别
    lm = fitlm(tbl, 'CBF~age+sex+group');

    % 提取统计数据
    mytstat(region) = lm.Coefficients{end, 'tStat'};
    mypval(region) = lm.Coefficients{end, 'pValue'};
end

sigregs=find(mypval<0.05);
pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values 如何计算出总体的P值？
sigregs_c = find(pvalue_fdr<0.05); % list 


% 显示结果

% 提取区域名称
regionNames = mergedData.Properties.VariableNames(2:380);
regionTypes = data4.regionType(1:379);
Networks = data4.network(1:379);
% 创建包含结果的表格
results_case_control = table(regionNames', mytstat, mypval, pvalue_fdr, regionTypes, Networks, 'VariableNames', {'RegionName', 'TStatistic', 'PValue', 'AdjustedPValue', 'RegionTypes', 'Networks'});

% 标记经过FDR矫正后显著的区域
results_case_control.SignificantAfterFDR = pvalue_fdr < 0.025;

% 筛选出显著的区域
significantRegions_case_control = results_case_control(results_case_control.SignificantAfterFDR, :);



cons = find(strcmp(group, 'HC'));

meancbf_con = mean(cbfData_raw(cons,:)',2);



pats = find(strcmp(group, 'CID'));

meancbf_case = mean(cbfData_raw(pats,:)',2);

% Calculate standard deviation for both groups

% 计算 meancbf_con 的平均值
averageValue_1 = mean(meancbf_con);

% 计算 meancbf_con 的标准差
stdDeviation_1 = std(meancbf_con);
% 计算 meancbf_con 的平均值
averageValue_2 = mean(meancbf_case);

% 计算 meancbf_con 的标准差
stdDeviation_2 = std(meancbf_case);

% Plot the correlation between the regional mean control MS and the mean t-statistic:为什么要做这一步？

%%进行spin test

addpath rotate_code/freesurfer_matalb/;
addpath rotate_code/rotate/;
addpath data;

% load('data/Shaefer.mat');
load('data/Glasser.mat');
% p_perm = perm_sphere_p(meancbf_con(1:360,:),mytstat(1:360,:),Glasser);
p_perm = perm_sphere_p(meancbf_case(1:360,:),mytstat(1:360,:),Glasser);
% 计算并显示相关系数
[r, p] = corr(meancbf_case, mytstat, 'Rows', 'complete');
corrText = sprintf('r = %.2f, $p_{\\mathrm{spin}} = %.3f$', r, p_perm);


% 绘制散点图
figure;
scatter(meancbf_case, mytstat, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7]);
% refline(0, 0); % 添加参考线
hold on;
coeffs = polyfit(meancbf_case, mytstat, 1);
fittedX = linspace(min(meancbf_case), max(meancbf_case), 200);
fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'Color', [1, 0.3, 0], 'LineWidth', 1.5);
hold off;
% plot([0 0], ylim, '-b'); % 添加垂直线
xlabel('Mean control CBF', 'FontName', 'Arial');
ylabel('Case-control t value', 'FontName', 'Arial');
% title('Correlation between control CBF value and t-statistics', 'FontName', 'Arial');
% legend(corrText, 'Location', 'best', 'FontName', 'Arial');
set(gca, 'FontSize', 20, 'FontName', 'Arial'); % 设置字体为 Arial
grid on; % 添加网格线


% 显示相关系数和p值
text(mean(meancbf_case), min(mytstat)-0.15, corrText, 'FontSize', 24, 'FontName', 'Arial', 'HorizontalAlignment', 'center','Interpreter', 'latex');
        % 保存图像
print('-dpng', fullfile('scatter_plots','Case_Control_t&control_corr.png'), '-r600');
close;
%%%计算残差和全局CBF值是否有差异

% 1. 计算残差
myresid_cbf = zeros(nregs, size(cbfData_raw, 1)); % 初始化残差矩阵
for region = 1:nregs
    y = cbfData_raw(:, region);
    [b, ~, resid] = regress(y, [ones(size(age)) age sex]);
    myresid_cbf(region, :) = resid;
end

% 2. 计算全局CBF均值
meanCBF = mean(myresid_cbf', 2); % 按行计算平均值

% 3. 评估组间差异
tbl_cbf = table(age, sex, group, meanCBF, 'VariableNames', {'Age', 'Sex', 'Group', 'MeanCBF'});
tbl_cbf.Sex = categorical(tbl_cbf.Sex);
tbl_cbf.Group = categorical(tbl_cbf.Group);
tbl_cbf.Group = categorical(tbl_cbf.Group, {'HC', 'CID'}); % 确保'HC'是第一个类别
lm_cbf = fitlm(tbl_cbf, 'MeanCBF~Age+Sex+Group');
p_global = lm_cbf.Coefficients{end, 'pValue'}; % 组别影响的p值

% 4. 可视化结果

pats = find(strcmp(group, 'CID'));

cons = find(strcmp(group, 'HC'));

x = reshape(myresid_cbf(:,cons),[nregs*length(cons),1]);
y = reshape(myresid_cbf(:,pats),[nregs*length(pats),1]);

% 残差图
figure;
h1 = histogram(x, 'FaceColor', [0 0.4470 0.7410]); % 深蓝色
hold on;
h2 = histogram(y, 'FaceColor', [0.8500 0.3250 0.0980]); % 橙色
h1.Normalization = 'probability';
h1.BinWidth = 2.00;
h2.Normalization = 'probability';
h2.BinWidth = 2.00;
legend('HCs', 'CID', 'Location', 'best', 'FontName', 'Arial');
xlim([-20.0 25.0]);
xlabel('Regional CBF residuals', 'FontName', 'Arial');
ylabel('Relative frequency', 'FontName', 'Arial');
set(gca, 'FontSize', 20, 'FontName', 'Arial'); % 设置字体为 Arial


% 保存图像
folder = 'scatter_plots'; % 指定文件夹
if ~exist(folder, 'dir')
    mkdir(folder); % 如果文件夹不存在，则创建
end

filename = fullfile(folder, 'frequency_residuals_plot.png'); % 指定文件名和路径
print('-dpng', '-r600', filename); % 保存为PNG格式，分辨率为300 DPI

% 全局CBF箱形图
figure;
b = boxplot(meanCBF, group, 'Labels', {'HC', 'CID'}, 'Notch', 'on', 'Whisker', 1.5);
title('Global mean CBF residuals by Group', 'FontName', 'Arial');
ylabel('Mean CBF residuals', 'FontName', 'Arial');
xlabel('Group', 'FontName', 'Arial');
set(gca, 'FontSize', 16, 'FontName', 'Arial'); % 设置字体为 Arial
set(findobj(gca, 'type', 'line'), 'Linewidth', 1.5); % 调整线宽


filename = fullfile(folder, 'meanCBF_residuals_plot.png'); % 指定文件名和路径
print('-dpng', '-r600', filename); % 保存为PNG格式，分辨率为300 DPI
% 输出组间差异的p值
disp(['P-value for the effect of group on global mean CBF: ', num2str(p_global)]);


%%%将结果映射到7-networks

networks=load('schaef7_overlap.txt'); % 
% 假设 networks 变量包含脑区分类信息
% 例如：networks = readmatrix('path_to_your_classification_data.csv');

% 初始化结果变量
mytstat_class = zeros(max(networks), 1);
mypval_class = zeros(max(networks), 1);

% 对每个脑区分类进行线性回归分析
for class = 1:max(networks)
    myregions = find(networks == class); % 找到属于当前分类的脑区
    classCBF = mean(cbfData_raw(:, myregions), 2); % 计算每个被试在该分类的平均CBF值

    % 创建表格并进行线性回归
    tbl = table(age, sex, group, classCBF, 'VariableNames', {'age', 'sex', 'group', 'ClassCBF'});
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    tbl.group = categorical(tbl.group, {'HC', 'CID'}); % 
    lm_class = fitlm(tbl, 'ClassCBF~age+sex+group');

    % 提取统计数据
    mytstat_class(class) = lm_class.Coefficients{end, 'tStat'};
    mypval_class(class) = lm_class.Coefficients{end, 'pValue'};
end

myp_class_fdr = mafdr(mypval_class,'BHFDR',1); % 
% 显示结果
network_results = table((1:max(networks))', mytstat_class, mypval_class, myp_class_fdr, 'VariableNames', {'Class', 'TStatistic', 'PValue', 'P_FDR'});




%%% 提取具有显著性区域的CBF值和量表数据

% 提取 CID 组的数据
cidData = mergedData(strcmp(mergedData.Group, 'CID'), :);
hcData = mergedData(strcmp(mergedData.Group, 'HC'), :);

% 确定显著区域的索引
significantRegions = significantRegions_case_control.RegionName;
sigRegionIndices = ismember(mergedData.Properties.VariableNames, significantRegions);

% 提取显著区域的 CBF 值
sigCBFData = cidData{:, sigRegionIndices};

% 为显著区域的 CBF 值设置正确的列名
sigCBFTable = array2table(sigCBFData, 'VariableNames', significantRegions);

% 提取 PSQI 和 ISI 值
psqiData = cidData.PSQI;
isiData = cidData.ISI;

% 合并所有数据
finalData = [table(cidData.SubjectID), sigCBFTable, table(psqiData, isiData)];

% 设置最终表格的列名
finalData.Properties.VariableNames = ['SubjectID', significantRegions', 'PSQI', 'ISI'];

% 指定输出文件名和路径
outputFilePath = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\code\sig_region&scales_data.xlsx';

% 保存数据到 Excel 文件
writetable(finalData, outputFilePath);

%%% 提取两个组的区域间CBF平均值
% 计算每个组在所有区域的平均值
% 排除非数值列
cidNumericData = cidData(:, varfun(@isnumeric, cidData, 'OutputFormat', 'uniform'));
hcNumericData = hcData(:, varfun(@isnumeric, hcData, 'OutputFormat', 'uniform'));

% 计算 CID 组各个区域的平均值
cidMeanValues = varfun(@mean, cidNumericData, 'OutputFormat', 'table');

% 计算 HC 组各个区域的平均值
hcMeanValues = varfun(@mean, hcNumericData, 'OutputFormat', 'table');

% 您可能还想给结果添加相应的标签
cidMeanValues.Group = repmat({'CID'}, size(cidMeanValues, 1), 1);
hcMeanValues.Group = repmat({'HC'}, size(hcMeanValues, 1), 1);

% 合并两组数据的平均值
groupMeans = [cidMeanValues; hcMeanValues];

% 指定输出文件名和路径
outputFilePathMeans = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\code\group_means_data.xlsx';

% 保存数据到 Excel 文件
writetable(groupMeans, outputFilePathMeans);





