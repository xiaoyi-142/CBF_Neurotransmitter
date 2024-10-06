clc; clear;

% 初始化数据
% Define the given values
n1 = 28; % Sample size of true acupuncture group
n2 = 28; % Sample size of sham acupuncture group
M1 = 5.93; % Mean change in PSQI score for true acupuncture group
M2 = 0.39; % Mean change in PSQI score for sham acupuncture group
SD1 = 2.36; % Standard deviation for true acupuncture group
SD2 = 2.02; % Standard deviation for sham acupuncture group
alpha = 0.05; % Two-tailed significance level
power = 0.95; % Desired statistical power
attrition_rate = 0.20; % Anticipated attrition rate

% Calculate the effect size 'd'
d = (M1 - M2) / sqrt(((n1 - 1) * SD1^2 + (n2 - 1) * SD2^2) / (n1 + n2 - 2));

% Z-scores for alpha and power
Z_alpha_2 = norminv(1 - alpha / 2);
Z_beta = norminv(power);

% Sample size calculation for each group
n = ((Z_alpha_2 + Z_beta)^2 * (SD1^2 + SD2^2)) / d^2;

% Adjusting for attrition rate
n_adjusted = n / (1 - attrition_rate);

% Rounding up to nearest whole number as you can't have a fraction of a subject
n_adjusted2 = ceil(n_adjusted);


% Display the result
fprintf('The adjusted sample size per group is: %d\n', n_adjusted2);














clc; clear;
% 设定参数
alpha = 0.05; % 显著性水平
beta = 0.20;  % 功效为80%，所以 beta = 1 - 0.80
effect_size = 1.0; % 假设的Cohen's d效应量

% 计算Z值
Z_alpha = norminv(1 - alpha / 2);  % 双尾测试的alpha
Z_beta = norminv(1 - beta);

% 样本量计算公式
n = ((Z_alpha + Z_beta)^2 * 2) / effect_size^2;

% 输出结果
fprintf('每组所需的样本量: %f\n', ceil(n)); % 向上取整








