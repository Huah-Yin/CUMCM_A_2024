clc;
clear;
close all;
params.mf = 4866; %浮子质量:kg
params.rf = 1; %浮子底半径：m
params.mz = 2433; %振子质量 (kg)
params.rou = 1025; %海水的密度 (kg/m3)
params.g = 9.8;
params.k4 = 656.3616; %垂荡兴波阻尼系数 (N·s/m)
params.m0 = 1335.535; %垂荡附加质量
params.k2 = 80000; %弹簧劲度系数
params.A = 6250; %垂荡激励力振幅 (N)
params.omega = 1.4005; %入射波浪频率 (s-1)

%%%%%%%线性方程组
params.detK = 10000; %直线阻尼器系数
%%%求解
dt = 0.001; %步长
t = 0:dt:179.4;
x0 = [0 0 0 0];
f = @(t, x) que1_1(t, x, params);

%%求解
[T, X] = ode45(f, t, x0);

%%竖直位移图
figure
plot(T, X(:, 1), 'r', 'LineWidth', 1.2); %浮子
hold on;
plot(T, X(:, 3), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); %振子
legend('浮子', '振子');
title("线性微分方程组 位移分析图");
xlabel('时间 (s)');
ylabel('位移 (m)');
grid on;

figure
plot(T, X(:, 1), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); %浮子
hold on;
plot(T, X(:, 1) - X(:, 3), 'r', 'LineWidth', 0.6); %差值
title('线性微分方程组 - 位移差值图');
legend('浮子', '振子');
xlabel('时间 (s)');
ylabel('位移 (m)');
grid on;

%%竖直速度图
figure;
plot(T, X(:, 2), 'r', 'LineWidth', 1.2); % 浮子，红色实线
hold on;
plot(T, X(:, 4), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); % 振子，蓝色虚线
legend('浮子', '振子');
title('线性微分方程组 - 速度分析图');
xlabel('时间 (s)');
ylabel('速度 (m/s)');
grid on;
%%浮子与振子速度差值图
figure;
plot(T, X(:, 2), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); % 浮子，红色实线
hold on;
plot(T, X(:, 2) - X(:, 4), 'r', 'LineWidth', 0.8); %差值
legend('浮子', '振子');
title("微分方程组为线性情况 速度差值分析图");
xlabel('时间 (s)');
ylabel('速度差值 (m/s)');
grid on;

% 获取这些时间点对应的位移和速度
y1 = X(:, 1); % 浮子位移
y2 = X(:, 2); % 浮子速度
y3 = X(:, 3); % 振子位移
y4 = X(:, 4); % 振子速度
% 计算位移差值和速度差值
diff_x = y1 - y3; % 位移差值
diff_v = y2 - y4; % 速度差值
% 准备写入 Excel 的数据
header = {'时间点 (s)', '浮子位移 (m)', '振子位移 (m)', '位移差值 (m)', '浮子速度 (m/s)', '振子速度 (m/s)', '速度差值 (m/s)'};
data = [T, y1, y3, diff_x, y2, y4, diff_v];

% 创建表格
Table = array2table(data, 'VariableNames', header);
% 将数据写入 Excel 文件
filename = 'F:\EXCEL\2022_a_1_1.xlsx';

writetable(Table, filename);
disp('数据已成功写入 2022_a_1_2.xlsx 文件');
