clear; clc; close all;
% 指定Excel文件的路径
filename = 'C:\Users\FEARLESS\Desktop\CID_ASL\source_data\post_data.xlsx';

% 读取工作表Sheet2
opts = detectImportOptions(filename, 'Sheet', 2);
data = readtable(filename, opts);

% 筛选RA和SA组
raData = data(strcmp(data.group, 'RA'), :);
saData = data(strcmp(data.group, 'SA'), :);

% 初始化结果数组
raResult = {};
saResult = {};

% 处理RA组数据
for i = 1:length(raData.Rawid)
    subjectData = raData(strcmp(raData.Rawid, raData.Rawid{i}), :);
    if sum(strcmp(subjectData.cond, 'Pre')) > 0 && sum(strcmp(subjectData.cond, 'Pos')) > 0
        raResult{end+1} = strjoin(unique(subjectData.SubjectID), ', ');
    end
end

% 处理SA组数据
for i = 1:length(saData.Rawid)
    subjectData = saData(strcmp(saData.Rawid, saData.Rawid{i}), :);
    if sum(strcmp(subjectData.cond, 'Pre')) > 0 && sum(strcmp(subjectData.cond, 'Pos')) > 0
        saResult{end+1} = strjoin(unique(subjectData.SubjectID), ', ');
    end
end
% 将结果转换为表格并添加编号
raTable = table((1:length(raResult))', raResult', 'VariableNames', {'Number', 'SubjectID'});
saTable = table((1:length(saResult))', saResult', 'VariableNames', {'Number', 'SubjectID'});

% 显示结果
disp('RA组治疗前后都有数据的被试ID:');
disp(raTable);

disp('SA组治疗前后都有数据的被试ID:');
disp(saTable);
