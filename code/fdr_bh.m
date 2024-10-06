% FDR校正函数
function adjusted_p = fdr_bh(pvalues, alpha)
    if nargin < 2
        alpha = 0.05; % 默认显著性水平
    end
    
    % 排序p值
    [sorted_p, order] = sort(pvalues);
    m = length(pvalues); % 总测试数
    
    % 计算调整后的p值
    adjusted_p = zeros(m, 1);
    prev_p = 0;
    for i = m:-1:1
        current_p = min(alpha * i / m, sorted_p(i));
        adjusted_p(order(i)) = max(current_p, prev_p);
        prev_p = adjusted_p(order(i));
    end
end