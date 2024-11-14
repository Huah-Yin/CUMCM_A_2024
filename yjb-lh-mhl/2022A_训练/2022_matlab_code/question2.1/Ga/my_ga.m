clc; clear; close all;
f1 = @(x) Max(x);
down = 37200; %下界
up =38000; %上界



[x, fval] = ga(f1, 1, [], [], [], [], down, up);
% 显示结果
disp(['最优解：', num2str(x)]);
disp(['最优目标值：', num2str(fval)]);
