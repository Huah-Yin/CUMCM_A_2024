function theta_result = delta(L, s, n)
    %用于求解两端点径向方向间初始夹角
    %输入:L 杆长
    %s 螺距宽
    %输出:两端点径向方向间初始夹角
    r = s * n;
    %     syms theta
    % 设置 fsolve 的选项
    options = optimoptions('fsolve', ...
        'TolX', 1e-12, ...                      % 变量变化精度
        'TolFun', 1e-12, ...   
         'Display', 'final', ...    % 函数值变化精度
        'MaxIterations', 1000, ...              % 最大迭代次数
        'MaxFunctionEvaluations', 5000, ...     % 最大函数评估次数
        'Algorithm', 'levenberg-marquardt');
    f = @(theta) cos(theta) - ((s * (2 * n * pi - theta) / (2 * pi)) ^ 2 + r ^ 2 - L ^ 2) / (2 * (s * (2 * n * pi - theta) / (2 * pi)) * r);

    % 使用 vpasolve 求解
    theta_result = fsolve(f, 0.5);
end
