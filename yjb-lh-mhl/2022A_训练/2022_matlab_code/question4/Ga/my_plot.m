clc; clear; close all;

tic;
% 设置范围和步长
x1 = 0:1000:100000; % X1 值
x2 = 0:1000:100000;

[x1, x2] = meshgrid(x1, x2);
% 计算每个网格点对应的函数值
PTO = arrayfun(@(a, b) Max([a, b]), x1, x2);

% 绘制3D网状图
mesh(x1, x2, PTO);

% 设置颜色映射和光照
colormap jet; % 使用 jet 颜色映射，色彩更丰富
colorbar; % 添加颜色条

% 调整视角
view(-30, 30); % 设置视角角度

% 添加轴标签和标题
xlabel('直线阻尼器系数', 'FontName', '宋体');
ylabel('旋转阻尼器系数', 'FontName', '宋体');
zlabel('PTO功率值', 'FontName', '宋体');
title('Max_p函数的3D网状图', 'FontName', '宋体');

toc;
