clc; clear; close all;

tic;
% 设置范围和步长
x1 = 1:1000:100000; % X1 值

% 计算每个 X1 对应的 Max 值
PTO = arrayfun(@(x) Max(x), x1);

% 绘制2D曲线图
figure;
plot(x1, PTO, 'LineWidth', 2); % 绘制曲线图并设置线宽

% 添加轴标签和标题
xlabel('X1 值');
ylabel('PTO 值');
title('Max 函数的2D曲线图');

% 添加网格线
grid on;

% 适应图形窗口
axis tight;

toc;