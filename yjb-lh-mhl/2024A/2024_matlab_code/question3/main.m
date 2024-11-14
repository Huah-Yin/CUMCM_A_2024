clc, clear, close all;
R = 4.5; % 最后的截止半径
b_max = 0.55;
b_min = 0.3;
dt = 1;
i = 1;
%%开始计时
tic;
while abs(b_max -b_min) > 0.01
    middle = (b_max + b_min) / 2;
    dx = middle; %螺距
    a = dx / (2 * pi);
    theta = R * 2 * pi / dx;
    n = floor(8.8 / dx) + 1;
    y0 = (a / 2) * ((2 * n * pi) * sqrt((2 * n * pi) ^ 2 + 1) + asinh(2 * n * pi)); %32 pi
    %%%%%得到时间末尾t
    t_ = y0 - (a / 2) * (theta * sqrt(theta ^ 2 + 1) + asinh(theta));
    t_end = floor(t_);
    my_flag = check(t_end, dx, dt, n); % % %返回1  没有发生碰撞

    if my_flag == 1 % %中间 b符合情况
        b_max = middle;
    else
        b_min = middle;
    end

    % max_ = zeros(1, 100);
    % min_ = zeros(1, 100);
    %     max_(i) = b_max;
    %     min_(i) = b_min;
    i = i + 1;
end
% 结束计时并获取时间
elapsedTime = toc;
disp(['二分法运行时间: ' num2str(elapsedTime) ' s']);
disp(i);
disp(dx);
