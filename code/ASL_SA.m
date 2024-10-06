% 清除环境
clear; clc; close all;

% 加载数据
filename = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\post_data.xlsx';
cbf_data = readtable(filename, 'Sheet', 'Post_CID_scheaf100');
ra_data = readtable(filename, 'Sheet', 'SA');

% 提取治疗前后的被试 SubjectID
pre_subjects = ra_data.SubjectID(strcmp(ra_data.Cond, 'Pre'));
pos_subjects = ra_data.SubjectID(strcmp(ra_data.Cond, 'Pos'));

% 提取脑区名称
region_names = cbf_data.Properties.VariableNames(2:end);
m = length(region_names); % 脑区数量

% 初始化存储结果的数组
t_values = zeros(m, 1);
p_values = zeros(m, 1);
t_values_significant = zeros(m, 1);
mean_pre = zeros(m, 1);  % 存储治疗前的平均值
mean_post = zeros(m, 1); % 存储治疗后的平均值

% 针对每个脑区执行配对 T 检验
for region = 1:m
    region_name = region_names{region};
    
    % 提取治疗前后的数据
    data_before = cbf_data{ismember(cbf_data.SubjectID, pre_subjects), region_name};
    data_after = cbf_data{ismember(cbf_data.SubjectID, pos_subjects), region_name};

    % 计算平均值
    mean_pre(region) = mean(data_before);
    mean_post(region) = mean(data_after);

    % 执行配对 T 检验
    [h, p, ci, stats] = ttest(data_before, data_after);
    
    % 保存结果
    t_values(region) = stats.tstat;
    p_values(region) = p;
end

% 对 p 值进行 FDR 校正
adjusted_p_values = mafdr(p_values, 'BHFDR', true);

% 将显著性结果保存到 t_values_significant
t_values_significant(p_values < 0.05) = t_values(p_values < 0.05);

% 保存结果为表格
result_table = table(region_names', t_values, p_values, adjusted_p_values, t_values_significant, mean_pre, mean_post, ...
                     'VariableNames', {'Region', 'TValue', 'PValue', 'AdjustedPValue', 'SigTValue', 'MeanPre', 'MeanPost'});

% 显示结果
disp(result_table);

% 将结果保存为Excel文件
output_filename = 'RCT_results\SA_result_table.xlsx';
writetable(result_table, output_filename);
