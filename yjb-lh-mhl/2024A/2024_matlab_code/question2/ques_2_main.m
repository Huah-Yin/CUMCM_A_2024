clc, clear, close all;
flag = 1;
z_k = 0.025; %板凳面积的占比系数
z = z_k * 0.66; %误差控制系数  m²   板凳面积的 3 %
number_i = 31; %从龙头紧接的板子（即就是 倒数第二块板子） 开始往后遍历个数

for t = 410:0.1:420

    %准备工作
    [r, theta, ~, ~] = prepare_2(t);
    %%%%%% alpha1 当作前把手时 ： 板方向与极矢方向的夹角
    %%%%%% alpha2 当作后把手时 ： 板方向与极矢方向的夹角
    [alpha1, alpha2] = Palpha(r);
    %特定角度
    beta0 = atan(15/27.5); % %板凳上任意一把手 坐标 与 右 or 左 顶点所连的直线 与 板凳方向的夹角

    if (beta0 > (pi / 2))
        exit;
    end

    %%求解龙头四顶点坐标式(x1,y1),(x3,y3),  对称后的笛卡尔坐标：(x2,y2),(x4,y4)
    r_1 = sqrt(r(1) ^ 2 + 0.15 ^ 2 + 0.275 ^ 2 - 2 * r(1) * sqrt(0.15 ^ 2 + 0.275 ^ 2) * cos(pi - alpha1(1) - beta0));
    %%%增加到 角度
    theta_1 = asin(sqrt(0.15 ^ 2 + 0.275 ^ 2) * sin(pi - alpha1(1) - beta0) / r_1);
    [x_1, y_1] = pol2cart(theta(1) - theta_1, r_1);

    [front_x, front_y] = pol2cart(theta(1), r(1));
    [behind_x, behind_y] = pol2cart(theta(2), r(2));

    [x_2, y_2] = symmetric_point([front_x, front_y], [behind_x, behind_y], [x_1, y_1]); % % %正确--x(1)

    r_3 = sqrt(r(2) ^ 2 + 0.15 ^ 2 + 0.275 ^ 2 - 2 * r(2) * sqrt(0.15 ^ 2 + 0.275 ^ 2) * cos(pi - alpha2(1) - beta0));
    theta_3 = asin(sqrt(0.15 ^ 2 + 0.275 ^ 2) * sin(pi - alpha2(1) - beta0) / r_3);
    [x_3, y_3] = pol2cart(theta(2) +theta_3, r_3); % % %正确

    [x_4, y_4] = symmetric_point([front_x, front_y], [behind_x, behind_y], [x_3, y_3]);

    %%求解 任意龙身四个顶点坐标(x1i,y1i),(x2i,y2i),(x3i,y3i),(x4i,y4i)
    x1 = zeros(30, 4);
    y1 = zeros(30, 4); % % %一个板凳四个坐标
    [L, ~] = size(x1); % % % % % TEST OK!

    for i = 2:number_i
        r1 = sqrt(r(i) ^ 2 + 0.15 ^ 2 + 0.275 ^ 2 - 2 * r(i) * sqrt(0.15 ^ 2 + 0.275 ^ 2) * cos(pi - alpha1(i) - beta0));
        theta1 = asin(sqrt(0.15 ^ 2 + 0.275 ^ 2) * sin(pi - alpha1(i) - beta0) / r1);
        [x1(i - 1, 1), y1(i - 1, 1)] = pol2cart(theta(i) - theta1, r1);

        [front_x2, front_y2] = pol2cart(theta(i), r(i));
        [behind_x2, behind_y2] = pol2cart(theta(i + 1), r(i + 1));
        % 通过1对称2
        [x1(i - 1, 2), y1(i - 1, 2)] = symmetric_point([front_x2, front_y2], [behind_x2, behind_y2], [x1(i - 1, 1), y1(i - 1, 1)]);

        r3 = sqrt(r(i + 1) ^ 2 + 0.15 ^ 2 + 0.275 ^ 2 - 2 * r(i + 1) * sqrt(0.15 ^ 2 + 0.275 ^ 2) * cos(pi - alpha2(i) - beta0));
        theta3 = asin(sqrt(0.15 ^ 2 + 0.275 ^ 2) * sin(pi - alpha2(i) - beta0) / r3);
        [x1(i - 1, 3), y1(i - 1, 3)] = pol2cart(theta(i + 1) + theta3, r3);
        %%%%% 通过3对称4
        [x1(i - 1, 4), y1(i - 1, 4)] = symmetric_point([front_x2, front_y2], [behind_x2, behind_y2], [x1(i - 1, 3), y1(i - 1, 3)]);
    end

    %%碰撞检测
    count_i = 2; % count_i为计数器:计数第几块板发生碰撞。

    for j = 2:L % % 因为：龙头始终会与其相邻的（即就是 倒数第二块板）相接，因此从倒数第三块板开始计数。
        %%%%%%%%%%% 龙头4个坐标与龙身坐标检测-面积法
        s1 = determin_picture(x_1, y_1, x1(j, :), y1(j, :));
        s2 = determin_picture(x_2, y_2, x1(j, :), y1(j, :));
        % s3 = determin_picture(x_3, y_3, x1(j, :), y1(j, :));
        % s4 = determin_picture(x_4, y_4, x1(j, :), y1(j, :));
        disp([t, s1, s2])

        if s1 <= z || s2 <= z % %|| s3 <= z || s4 <= z
            flag = 0;
            disp([count_i + 2, s1, s2]);
            % disp([count_i, s1, s2, s3, s4]);
            break;
        end

        count_i = count_i + 1;
    end

    %信号判断
    if flag == 0 % % % % %已经发生了碰撞，终止循环，得出此时此刻的所有v 坐标
        [~, ~, ~, v] = prepare_2(t - 1);
        v = v';
        disp("在ts 发生了碰撞");
        disp(t);
        break;
    end

end

if flag == 1
    disp('增加时间');
end

[x, y] = pol2cart(theta, r);
