% 清除环境
clear; clc; close all;

% 加载数据
filename = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\post_data.xlsx';
cbf_data = readtable(filename, 'Sheet', 'Post_CID_scheaf100');
ra_data = readtable(filename, 'Sheet', 'RA');

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
output_filename = 'RCT_results\RA_result_table.xlsx';
writetable(result_table, output_filename);


%%进行相关性评估

% 找出显著性脑区
significant_regions = region_names(p_values < 0.05);

% 初始化存储相关性分析结果的结构
correlation_results = struct;

% 遍历每个显著的脑区
for i = 1:length(region_names)

    % 提取该脑区的治疗前后数据
    region = region_names{i};
    data_before = cbf_data{ismember(cbf_data.SubjectID, pre_subjects), region};
    data_after = cbf_data{ismember(cbf_data.SubjectID, pos_subjects), region};

    % 计算CBF改变值
    cbf_change = data_after - data_before;

    % 提取临床量表的改变值
    clinical_metrics = {'PSQI', 'ISI', 'SAS', 'SDS'};
    for j = 1:length(clinical_metrics)
        metric = clinical_metrics{j};
        pre_metric = ra_data{strcmp(ra_data.Cond, 'Pre') & ismember(ra_data.SubjectID, pre_subjects), metric};
        post_metric = ra_data{strcmp(ra_data.Cond, 'Pos') & ismember(ra_data.SubjectID, pos_subjects), metric};

        clinical_change = post_metric - pre_metric;

        % 相关性分析
        [R, P] = corrcoef(cbf_change, clinical_change, 'Rows', 'complete');

        % 保存结果
        correlation_results.(region).(metric).R = R(1,2);
        correlation_results.(region).(metric).P = P(1,2);
    end
end


% 显示相关性分析结果
disp(correlation_results);



% 定义显著性水平
significance_level = 0.05;

% 初始化用于保存显著相关结果的结构体
significant_correlations = struct;

% 遍历 correlation_results 结构体
regions = fieldnames(correlation_results);
for i = 1:length(regions)
    region = regions{i};
    metrics = fieldnames(correlation_results.(region));
    for j = 1:length(metrics)
        metric = metrics{j};
        if correlation_results.(region).(metric).P < significance_level
            % 如果 p 值小于显著性水平，保存结果
            significant_correlations.(region).(metric) = correlation_results.(region).(metric);
        end
    end
end

% 显示具有显著相关性的结果
disp(significant_correlations);





% 定义显著性水平
significance_level = 0.05;

% 遍历每个脑区
for i = 1:length(region_names)
    region = region_names{i};

    % 遍历每个临床量表
    for j = 1:length(clinical_metrics)
        metric = clinical_metrics{j};

        % 检查相关性是否显著
        if correlation_results.(region).(metric).P < significance_level
            % 提取相关数据
            data_before = cbf_data{ismember(cbf_data.SubjectID, pre_subjects), region};
            data_after = cbf_data{ismember(cbf_data.SubjectID, pos_subjects), region};
            cbf_change = data_after - data_before;
            
            pre_metric = ra_data{strcmp(ra_data.Cond, 'Pre') & ismember(ra_data.SubjectID, pre_subjects), metric};
            post_metric = ra_data{strcmp(ra_data.Cond, 'Pos') & ismember(ra_data.SubjectID, pos_subjects), metric};
            clinical_change = post_metric - pre_metric;

            % 绘制散点图
            figure;
            scatter(clinical_change, cbf_change, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0.7 0.7 0.7]);
            
            % 回归线
            hold on;
            coeffs = polyfit(clinical_change, cbf_change, 1);
            fittedX = linspace(min(clinical_change), max(clinical_change), 200);
            fittedY = polyval(coeffs, fittedX);
            plot(fittedX, fittedY, 'Color', [1, 0.3, 0], 'LineWidth', 1.5);
            hold off;

            % 标注相关性和 P 值
            text_pos_x = min(clinical_change) + 0.25 * range(clinical_change);
            text_pos_y = min(cbf_change) + 0.10 * range(cbf_change);
            text(text_pos_x, text_pos_y, sprintf('r = %.2f, p = %.3f', correlation_results.(region).(metric).R, correlation_results.(region).(metric).P), 'FontSize', 18);
            
            region = strrep(region, '_', ' ');
            % 标签和标题
            xlabel([metric ' changes'], 'FontSize', 18);
            ylabel([region ' CBF changes'], 'FontSize', 18);
%             title(['Correlation between ' metric ' and ' region ' CBF changes']);
            set(gca, 'FontSize', 18, 'FontName', 'Arial'); 
            % 保存图像
            saveas(gcf, fullfile('RCT_results', [metric '_', strrep(region, ' ', '_'), '.png']));
            close;
        end
    end
end
