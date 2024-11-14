clc;
clear;
close all;
params.mf = 4866; %浮子质量:kg
params.rf = 1; %浮子底半径：m
params.mz = 2433; %振子质量 (kg)
params.rou = 1025; %海水的密度 (kg/m3)
params.g = 9.8;
params.k4 = 528.5018; %垂荡兴波阻尼系数 (N·s/m)
params.m0 = 1091.099; %垂荡附加质量
params.k2 = 80000; %弹簧劲度系数
params.A = 1760; %垂荡激励力振幅 (N)
params.omega = 1.9806; %入射波浪频率 (s-1)

params.L1 = 8890.7;
params.L2 = 250000;
params.L3 = 1000; %和旋转阻尼器的阻尼系数
params.L4 = 1655.909;
params.L = 2140;
params.Jc = 7142.493;
params.o = 0.764;
%%%%%%%线性方程组
params.detK = 10000; %直线阻尼器系数
%%%求解
dt = 0.2; %步长
t = 0:dt:179.4;
x0 = [0 0 0 0 0 0 0 0];
f = @(t, x) que3(t, x, params);

%%求解
[T, X] = ode45(f, t, x0);

%竖直位移图
figure
plot(T, X(:, 1), 'r', 'LineWidth', 1.2); %浮子
hold on;
plot(T, X(:, 3), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); %振子
legend('浮子', '振子', 'FontName', '宋体');
title("线性微分方程组 位移分析图", 'FontName', '宋体');
xlabel('时间 (s)', 'FontName', '宋体');
ylabel('位移 (m)', 'FontName', '宋体');
grid on;

figure
plot(T, X(:, 1), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); %浮子
hold on;
plot(T, X(:, 1) - X(:, 3), 'r', 'LineWidth', 0.6); %差值
title('线性微分方程组 - 位移差值图', 'FontName', '宋体');
legend('浮子', '振子', 'FontName', '宋体');
xlabel('时间 (s)', 'FontName', '宋体');
ylabel('位移 (m)', 'FontName', '宋体');
grid on;

%%竖直速度图
figure;
plot(T, X(:, 2), 'r', 'LineWidth', 1.2); % 浮子，红色实线
hold on;
plot(T, X(:, 4), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); % 振子，蓝色虚线
legend('浮子', '振子', 'FontName', '宋体');
title('线性微分方程组 - 速度分析图', 'FontName', '宋体');
xlabel('时间 (s)', 'FontName', '宋体');
ylabel('速度 (m/s)', 'FontName', '宋体');
grid on;
%%浮子与振子速度差值图
figure;
plot(T, X(:, 2), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); % 浮子，红色实线
hold on;
plot(T, X(:, 2) - X(:, 4), 'r', 'LineWidth', 0.8); %差值
legend('浮子', '振子', 'FontName', '宋体');
title("微分方程组为线性情况 速度差值分析图", 'FontName', '宋体');
xlabel('时间 (s)', 'FontName', '宋体');
ylabel('速度差值 (m/s)', 'FontName', '宋体');
grid on;

%%%%%%%%%%纵摇
%角位移
figure
plot(T, X(:, 5), 'r', 'LineWidth', 1.2); %浮子
hold on;
plot(T, X(:, 7), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); %振子
legend('浮子', '振子', 'FontName', '宋体');
title("线性微分方程组角位移分析图", 'FontName', '宋体');
xlabel('时间/s', 'FontName', '宋体');
ylabel('角位移/rad', 'FontName', '宋体');
grid on;
%角位移之差
figure
plot(T, X(:, 5), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); %浮子
hold on;
plot(T, X(:, 5) - X(:, 7), 'r', 'LineWidth', 0.6); %差值
title('线性微分方程组角位移差值图', 'FontName', '宋体');
legend('浮子', '振子', 'FontName', '宋体');
xlabel('时间/s', 'FontName', '宋体');
ylabel('角位移/rad', 'FontName', '宋体');
grid on;

%%角速度
figure;
plot(T, X(:, 6), 'r', 'LineWidth', 1.2); % 浮子，红色实线
hold on;
plot(T, X(:, 8), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); % 振子，蓝色虚线
legend('浮子', '振子', 'FontName', '宋体');
title('线性微分方程组 角速度分析图', 'FontName', '宋体');
xlabel('时间 s', 'FontName', '宋体');
ylabel('角速度 rad/s', 'FontName', '宋体');
grid on;
%%浮子与振子角速度差值图
figure;
plot(T, X(:, 6), 'Color', [0.6 0.8 1], 'LineWidth', 0.8); % 浮子，红色实线
hold on;
plot(T, X(:, 6) - X(:, 8), 'r', 'LineWidth', 0.8); %差值
legend('浮子', '振子', 'FontName', '宋体');
title("微分方程组为线性情况 角速度差值分析图", 'FontName', '宋体');
xlabel('时间 s', 'FontName', '宋体');
ylabel('角速度差值 rad/s', 'FontName', '宋体');
grid on;

% 提取特定时间点的数据
time_points = [10, 20, 40, 60, 100];
num_points = length(time_points);

% 预分配数组
selected_data = zeros(num_points, 7);

for i = 1:num_points
    % 找到最接近的时间点的索引
    [~, idx] = min(abs(T - time_points(i)));

    % 提取数据
    selected_data(i, 1) = T(idx); % 时间点
    selected_data(i, 2) = X(idx, 1); % 浮子角位移
    selected_data(i, 3) = X(idx, 3); % 振子角位移
    selected_data(i, 4) = X(idx, 1) - X(idx, 3); % 角位移差值
    selected_data(i, 5) = X(idx, 2); % 浮子角速度
    selected_data(i, 6) = X(idx, 4); % 振子角速度
    selected_data(i, 7) = X(idx, 2) - X(idx, 4); % 角速度差值
end

% 创建表格
header = {'时间点 (s)', '浮子角位移 rad', '振子角位移 rad', '角位移差值 rad', '浮子角速度 rad/s', '振子角速度 rad/s', '角速度差值 rad/s'};
Table = array2table(selected_data, 'VariableNames', header);

% 将数据写入 Excel 文件
filename = 'F:\EXCEL\que_3_1.xlsx';
writetable(Table, filename);

disp('特定时间点的数据已成功写入 que_3_1.xlsx 文件');
