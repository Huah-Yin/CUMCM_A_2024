clc, clear, close all;
%%初始化参数
L1 = 2.86; %龙头长度 单位：m
L2 = 1.65; %龙身龙尾长度
s = 0.55; %初始螺线宽距
b = s / (2 * pi);
delta_theta_1 = delta(L1, s);
delta_theta_2 = delta(L2, s);
dt = 0.1;
t_end = 352;
t0 = 0:dt:t_end;
t0 = t0';
len_t = length(t0);
v = zeros(len_t, 224); %定义线速度,一共224个把手
theta = v;
Omega = v; %定义角速度
r = v; %定义极径
v(:, 1) = 1;
%%计算第一个点的角速度

w1_funs = @(t, theta) 1 ./ (b * sqrt(1 + (32 * pi - theta) .^ 2));
theta(:, 1) = 32 * pi - dis_euler(w1_funs, 1, t_end * 10, 0);
r(:, 1) = b .* theta(:, 1);
Omega(:, 1) = 1 ./ sqrt(r(:, 1) .^ 2 + b ^ 2);

%%计算第二个的角速度
i = 2 * ones(1, 222);

while (abs(theta(i(1), 1) - theta(1, 1)) < delta_theta_1)
    i (1) = i(1) + 1;
end

%刚性杆的约束条件
calpha1 = @(t, theta_2) (L1 ^ 2 + r(t, 1) .^ 2 - (b * (32 * pi - theta_2)) .^ 2) ./ (2 * L1 .* r(t, 1));
salpha1 = @(t, theta_2)sqrt(1 - calpha1(t, theta_2) .^ 2);
calpha2 = @(t, theta_2)((b * (32 * pi - theta_2)) .^ 2 + L1 ^ 2 - r(t, 1) .^ 2) ./ (2 * b .* (32 * pi - theta_2) .* L1);
salpha2 = @(t, theta_2)sqrt(1 - calpha2(t, theta_2) .^ 2);
% 沿杆的分速度相等
funs = @(t, theta_2) Omega(t, 1) .* (calpha1(t, theta_2) + theta(t, 1) .* salpha1(t, theta_2)) ./ (calpha2(t, theta_2) + (32 * pi - theta_2) .* salpha2(t, theta_2));
y0 = 0;
h = 0.1;
theta(i(1):end, 2) = 32 * pi - dis_euler(funs, i(1), 3000, y0);
r(i(1):end, 2) = b .* theta(i(1):end, 2);

for k = i(1):length(t0)
    Omega(k, 2) = funs(k, 32 * pi - theta(k, 2));
end

v(i(1):end, 2) = Omega(i(1):end, 2) .* sqrt(r(i(1):end, 2) .^ 2 + b ^ 2);
%计算第二个点以后的角速度
flag = 0;

for j = 3:224
    i(j - 1) = i(j - 2);

    while (abs(theta(i(j - 1) + 1, j - 1) - theta(i(j - 2) + 1, j - 1)) < delta_theta_2)
        i (j - 1) = i(j - 1) + 1;

        if i(j - 1) > 3000
            flag = 1;
            break;
        end

    end

    if flag == 1
        break;
    end

    calpha1 = @(t, theta_2) (L2 ^ 2 + r(t, j - 1) .^ 2 - (b * (32 * pi - theta_2)) .^ 2) ./ (2 * L2 .* r(t, j - 1));
    salpha1 = @(t, theta_2)sqrt(1 - calpha1(t, theta_2) .^ 2);
    calpha2 = @(t, theta_2)((b * (32 * pi - theta_2)) .^ 2 + L2 ^ 2 - r(t, j - 1) .^ 2) ./ (2 * b .* (32 * pi - theta_2) .* L2);
    salpha2 = @(t, theta_2)sqrt(1 - calpha2(t, theta_2) .^ 2);
    funs = @(t, theta_2) Omega(t, j - 1) .* (calpha1(t, theta_2) + theta(t, j - 1) .* salpha1(t, theta_2)) ./ (calpha2(t, theta_2) + (32 * pi - theta_2) .* salpha2(t, theta_2));
    y0 = 0;
    h = 0.1;
    theta(i(j - 1):end, j) = 32 * pi - dis_euler(funs, i(j - 1), 3000, y0);
    r(i(j - 1):end, j) = b .* theta(i(j - 1):end, j);

    for k = i(j - 1):length(t0)
        Omega(k, j) = funs(k, 32 * pi - theta(k, j));
    end

    v(i(j - 1):end, j) = Omega(i(j - 1):end, j) .* sqrt(r(i(j - 1):end, j) .^ 2 + b ^ 2);
end

%%转化为直角坐标
[x, y] = pol2cart(theta, r);
y(1, 1) = 0;

for n = 1:222

    if i(n) <= 3001
        y(i(n), n + 1) = 0;
    else
        break;
    end

end

%%%%将没进入盘龙的板凳坐标更新为 nan
for i = 2:224

    for j = 1:len_t - 1

        y(j, i) = nan;

        if y(j + 1, i) ~= 0
            y(j, i) = 0;
            break;
        else
            y(j + 1, i) = nan;
        end

    end

end

for i = 2:224

    for j = 1:len_t - 1

        x(j, i) = nan;

        if x(j + 1, i) == 0
            x(j, i) = nan;

            if j == len_t - 1
                x(j + 1, i) = nan;
            end

            continue;
        else
            x(j, i) = nan;
            break;
        end

    end

end

%%%% 螺旋线-三维空间
t_step = 0:1:300;
t_step = (t_step / dt) + 1;
figure;
plot3(t0(t_step), x(t_step), y(t_step), 'LineWidth', 2); % 使用 plot3 绘制三维曲线
grid on; % 打开网格
title('三维空间中的曲线 (x, y, t)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('时间 (t)', 'FontSize', 12, 'FontWeight', 'bold'); % x 轴为时间 t
ylabel('X 坐标', 'FontSize', 12, 'FontWeight', 'bold'); % y 轴为 x
zlabel('Y 坐标', 'FontSize', 12, 'FontWeight', 'bold'); % z 轴为 y
set(gca, 'FontSize', 12); % 设置坐标轴字体大小

%%%%%速度的三维空间

% figure;
% handles = 1:1:224; % 把手编号

% % 生成时间和把手编号的网格
% [H, T] = meshgrid(handles, t0(t_step));

% % 绘制网格图
% mesh(H, T, v(t_step, :));

% % 设置图形标题和轴标签
% title('把手速度随时间的变化', 'FontSize', 14, 'FontWeight', 'bold');
% xlabel('把手编号', 'FontSize', 12, 'FontWeight', 'bold'); % x 轴为把手编号
% ylabel('时间 (s)', 'FontSize', 12, 'FontWeight', 'bold'); % y 轴为时间
% zlabel('速度 (m/s)', 'FontSize', 12, 'FontWeight', 'bold'); % z 轴为速度

% % 打开网格，调整坐标轴
% grid on;
% set(gca, 'FontSize', 12); % 设置坐标轴字体大小

% % 去掉 shading 插值（网格不需要插值）
% colorbar; % 添加颜色条，表示速度大小

% 绘制三维曲面图
figure;
handles = 1:1:224;

[H, T] = meshgrid(handles, t0(t_step));
surf(H, T, v(t_step, :));
% 设置图形标题和轴标签
title('把手速度随时间的变化', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('把手编号', 'FontSize', 12, 'FontWeight', 'bold'); % x 轴为把手编号
ylabel('时间 (s)', 'FontSize', 12, 'FontWeight', 'bold'); % y 轴为时间
zlabel('速度 (m/s)', 'FontSize', 12, 'FontWeight', 'bold'); % z 轴为速度

% 打开网格，调整坐标轴
grid on;
set(gca, 'FontSize', 12); % 设置坐标轴字体大小
shading interp; % 平滑曲面插值，去除网格线
colorbar; % 添加颜色条，表示速度大小
% % 绘制第 300 秒的散点图
% figure;
% scatter(x(end, :), y(end, :), 50, 'filled'); % 绘制散点图，'50' 是点的大小，'filled' 表示填充点
% title('300 秒时的散点图', 'FontSize', 14, 'FontWeight', 'bold'); % 设置标题
% xlabel('X 坐标', 'FontSize', 12, 'FontWeight', 'bold'); % 设置 x 轴标签
% ylabel('Y 坐标', 'FontSize', 12, 'FontWeight', 'bold'); % 设置 y 轴标签
% grid on; % 打开网格
% set(gca, 'FontSize', 12); % 设置坐标轴字体大小
% axis equal; % 保持 x 和 y 轴比例相等
% hold on
% plot(x(end, :), y(end, :))

x_tran = x'; % 转置后的 x
y_tran = y'; % 转置后的 y
v_tran = v'; % 转置后的 v
step = 0:1:300;
step_index = step / dt +1;
% 对 x, y, v 进行采样
x_1 = x_tran(:, step_index); % 按时间间隔采样 x
y_1 = y_tran(:, step_index);
v_1 = v_tran(:, step_index);

excel_position = zeros(224 * 2, 301);
excel_v = zeros(224, 301);

for i = 1:224 * 2

    if mod(i, 2) ~= 0 %为奇数
        j = floor(i / 2) + 1; % (i/2 ) +1
        excel_position(i, :) = x_1(j, :);
    else %偶数
        k = i / 2;
        excel_position(i, :) = y_1(k, :);
    end

end

% %%%%%%%%%%%%论文当中的表
% colum_step = [1 61 121 181 241 301]; %下标
% row_step1 = [0 1 51 101 151 201 223] + 1;
% row_step2 = [1 2 3 4 103 104 203 204 303 304 403 404 447 448];
% paper_matrix1 = excel_position(row_step2, colum_step); %position
% paper_matrix2 = v_1(row_step1, colum_step); %speed
% header = {' 0 s', '60 s', ' 120 s', '180 s', '240 s', '300 s'};
% Table1 = array2table(paper_matrix1, 'VariableNames', header);
% Table2 = array2table(paper_matrix2, 'VariableNames', header);
% % 将数据写入 Excel 文件
% filename1 = 'F:\EXCEL\position.xlsx';
% filename2 = 'F:\EXCEL\speed.xlsx';
% writetable(Table1, filename1);
% writetable(Table2, filename2);
% disp('paper 速度 and 位置 已写入!');

% "F:\EXCEL\result1.xlsx"
%%%%%%%写入文件

% 文件路径
filename = 'F:\EXCEL\result1.xlsx';

[~, ~, raw_data] = xlsread(filename, '位置');
start_row = 2;
start_col = 2;

[row, col] = size(excel_position);
raw_data(start_row:(start_row + row - 1), start_col:(start_col + col - 1)) = num2cell(excel_position);

writecell(raw_data, filename, 'Sheet', '位置');

disp('矩阵已成功写入指定位置，并保留了文字信息');
writetable(excel_position, filename, 'Sheet', '位置');

disp('x, y 已写入 "位置" 表，速度 v 已写入 "速度" 表');
