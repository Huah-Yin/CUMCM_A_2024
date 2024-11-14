clc; clear; close all;
format long g
dt = 1;
dx = 1.7; %m
v0 = 1;
a = dx / (2 * pi);
R = (2.83 * 3) / 2;
n = round(8.8 / dx);
theta = R * 2 * pi / dx;
y0 = (a / 2) * ((2 * n * pi) * sqrt((2 * n * pi) ^ 2 + 1) + asinh(2 * n * pi));
t_ = y0 - (a / 2) * (theta * sqrt(theta ^ 2 + 1) + asinh(theta));
t_end = floor(t_);
[r0, theta0, Omega, v_0] = prepare_3(t_end, dx, dt, n);

[x, y] = pol2cart(theta0, r0);
% x = x';
% y = y';
% v_0 = v_0';
% scatter(x, y);

% %%%%%%%%%%%%论文当中的表
% step = [1 51 101 151 201] + 1; %下标
% paper_matrix1 = x(step); %#ok<> %x
% paper_matrix2 = y(step); %y
% paper_matrix3 = v_0(step);
% combined_matrix = [paper_matrix1, paper_matrix2, paper_matrix3];
% header = {'x', 'y', 'v'};
% Table1 = array2table(combined_matrix, 'VariableNames', header);
% % 将数据写入 Excel 文件
% filename1 = 'F:\EXCEL\que2_position_v.xlsx';
% % filename2 = 'F:\EXCEL\speed.xlsx';
% writetable(Table1, filename1);
% disp('paper 速度 and 位置 已写入!');

% step = 0:1:100;
% step_index = step / dt +1;
% % 对 x, y, v 进行采样
% x_1 = x(:, step_index); % 按时间间隔采样 x
% y_1 = y(:, step_index);
% v_1 = v_0(:, step_index);

x_1 = x';
y_1 = y';
v_1 = v_0';
excel_position = zeros(224 * 2, 101);
excel_v = zeros(224, 101);

for i = 1:224 * 2

    if mod(i, 2) ~= 0 %为奇数
        j = floor(i / 2) + 1; % (i/2 ) +1
        excel_position(i, :) = x_1(j, :);
    else %偶数
        k = i / 2;
        excel_position(i, :) = y_1(k, :);
    end

end

% 文件路径
filename = 'F:\EXCEL\result4.xlsx';

[~, ~, raw_data] = xlsread(filename, '位置');
start_row = 2;
start_col = 2;

[row, col] = size(excel_position);
raw_data(start_row:(start_row + row - 1), start_col:(start_col + col - 1)) = num2cell(excel_position);

writecell(raw_data, filename, 'Sheet', '位置');

disp('矩阵已成功写入指定位置，并保留了文字信息');
writetable(excel_position, filename, 'Sheet', '位置');

disp('x, y 已写入 "位置" 表，速度 v 已写入 "速度" 表');
