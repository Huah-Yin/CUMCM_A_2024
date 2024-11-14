clc; clear; close all;
%%本文件针对于第一个版本比较欧拉法与解析式法的误差矩阵
format long g
%%% 参数
dt = 0.1;
t = 0:dt:300;
len_t = length(t);
cir_num = 16; % 16 圈
dx = 0.55; % 间隔 55cm
len_A = cir_num * dx;

%%%%%%%%%    r = a * theta
a = dx / (2 * pi);

%%%%%%%%%
theta_result = zeros(1, len_t);
diff_ = zeros(1, len_t);

% 初始猜测值
theta_guess = 32 * pi;

% 初始值 y0 的计算
y0 = (a / 2) * ((32 * pi) * sqrt((32 * pi) ^ 2 + 1) + asinh(32 * pi));

% 设置 fsolve 的高精度选项
options = optimoptions('fsolve', ...
    'Display', 'off', ... % 关闭迭代信息显示
    'TolFun', 1e-7, ... % 函数值容差
    'TolX', 1e-7, ... % 解的容差
    'MaxIter', 1000, ... % 最大迭代次数
    'MaxFunEvals', 2000); % 最大函数计算次数

for i = 1:len_t
    % 定义目标函数作为 fsolve 输入
    fun = @(theta) y0 - (a / 2) * (theta * sqrt(theta ^ 2 + 1) + asinh(theta)) - t(i);

    % 使用 fsolve 求解 theta
    [theta_result(i), diff_(i)] = fsolve(fun, theta_guess, options);

    % 更新猜测值（步长调整为 -0.01）
    theta_guess = theta_result(i) - 0.01;
end

theta_colum = theta_result'; %转置矩阵
%%初始化参数
L1 = 2.86; %龙头长度 单位：m
L2 = 1.65; %龙身龙尾长度
s = 0.55; %初始螺线宽距
b = s / (2 * pi);
delta_theta_1 = delta(L1, s);
delta_theta_2 = delta(L2, s);
t0 = 0:0.1:300;
t0 = t0';
v = zeros(length(t0), 1); %定义线速度
theta = v;
w = v; %定义角速度
r = v; %定义极径
v(:) = 1;
%%计算第一个点的角速度

w1_funs = @(t, theta) 1 ./ (b * sqrt(1 + (32 * pi - theta) .^ 2));
theta = 32 * pi - dis_euler(w1_funs, 1, 3000, 0);
r = b .* theta;
w = 1 ./ sqrt(r .^ 2 + b ^ 2);

deltaz = abs(theta - theta_colum); %欧拉法和解析式法误差矩阵
element =  power((theta - theta_colum) , 2);
mse = mean(element);


% 绘制逐元素 MSE 图
figure;
plot(element,  'LineWidth', 2);
title(['逐元素均方误差 (MSE)，总 MSE = ', num2str(mse)]);
xlabel('元素索引');
ylabel('平方误差');
grid on;
disp("误差平方和是 : ");
figure;
boxplot(element);  % 绘制箱线图
title('逐元素误差平方的箱线图', 'FontSize', 14, 'FontWeight', 'bold');  % 设置标题
xlabel('数据组', 'FontSize', 12, 'FontWeight', 'bold');  % 设置 x 轴标签
ylabel('误差平方值', 'FontSize', 12, 'FontWeight', 'bold');  % 设置 y 轴标签
set(gca, 'FontSize', 12);  % 设置坐标轴字体大小

grid on;  % 打开网格
z = sum(element); %误差平方和
disp(z);
