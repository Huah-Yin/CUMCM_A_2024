clc; clear; close all;
f1 = @(x) Max(x);
down = [60000 , 80000]; %下界
up = [100000 ,100000]; %上界

[x, fval] = ga(f1, 2, [], [], [], [], down, up);
% 显示结果
disp(['最优解：', num2str(x)]);
disp(['最优目标值：', num2str(fval)]);
