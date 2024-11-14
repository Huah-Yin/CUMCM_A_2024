clc; clear; close all;
f1 = @(x) Max(x);
down = [80000, 0]; %下界
up = [100000, 0.5]; %上界

% 设置绘图选项，'gaplotbestf' 用于绘制迭代过程中最佳目标值的图
options = optimoptions('ga', 'PlotFcn', @gaplotbestf);

[x, fval] = ga(f1, 2, [], [], [], [], down, up,[],options);
% 显示结果
disp(['最优解：', num2str(x)]);
disp(['最优目标值：', num2str(fval)]);
