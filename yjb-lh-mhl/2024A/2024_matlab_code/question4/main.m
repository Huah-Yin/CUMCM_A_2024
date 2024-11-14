clc; clear; close all;
format long g
dt = 1;
dx = 1.7; %m
v0 = 1;
i = 1;
a = dx / (2 * pi);
r_max = 4.5;
r_min = (2.86 * 3) / 2; % % % % % %需对R进行特殊额条件判断
%%%%%%%二分搜索->找最小的r
while r_max - r_min > 0.00001
    middle = (r_max + r_min) / 2;
    R = middle; % % % % %将middle 赋值给中间区域

    theta = R * 2 * pi / dx;
    n = round(8.8 / dx);
    y0 = (a / 2) * ((2 * n * pi) * sqrt((2 * n * pi) ^ 2 + 1) + asinh(2 * n * pi)); %32 pi
    %%%%%得到时间末尾t
    t_ = y0 - (a / 2) * (theta * sqrt(theta ^ 2 + 1) + asinh(theta));
    t_end = floor(t_);
    my_flag = check(t_end, dx, dt, n); % % %返回1  没有发生碰撞

    if my_flag == 1 % %中间 r符合情况
        r_max = middle;
    else
        r_min = middle;
    end

    i = i + 1;
end

disp(i);
disp(R);
disp(r_min);
disp(r_min);
