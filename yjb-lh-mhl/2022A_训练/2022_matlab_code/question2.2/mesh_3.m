clc; clear; close all;

tic;
% 设置范围和步长
x1 = 0:10000:100000; % 第一列的值
x2 = 0:0.01:1; % 第二列的值

[x1, x2] = meshgrid(x1, x2);
% 计算每个网格点对应的函数值
Max_p = arrayfun(@(a, b) Max([a, b]), x1, x2);

% 绘制3D网状图
mesh(x1, x2, Max_p);

% 设置颜色映射和光照
colormap jet; % 使用 jet 颜色映射，色彩更丰富
colorbar; % 添加颜色条

% 调整视角
view(-30, 30); % 设置视角角度

% 添加轴标签和标题
xlabel('X1 值');
ylabel('X2 值');
zlabel('Max 值');
title('Max 函数的3D网状图');

toc;


