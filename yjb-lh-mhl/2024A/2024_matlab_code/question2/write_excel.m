clc; clear; close all;

[r0, theta0, Omega, v0] = prepare_2(412);

[x, y] = pol2cart(theta0, r0);
x = x';
y = y';
v0 = v0';
scatter(x, y);
% %%%%%%%%%%%%论文当中的表
step = [1 51 101 151 201] + 1; %下标
paper_matrix1 = x(step); %x
paper_matrix2 = y(step); %y
paper_matrix3 = v0(step);
combined_matrix = [paper_matrix1, paper_matrix2, paper_matrix3];
header = {'x', 'y', 'v'};
Table1 = array2table(combined_matrix, 'VariableNames', header);
% 将数据写入 Excel 文件
filename1 = 'F:\EXCEL\que2_position_v.xlsx';
% filename2 = 'F:\EXCEL\speed.xlsx';
writetable(Table1, filename1);
disp('paper 速度 and 位置 已写入!');
