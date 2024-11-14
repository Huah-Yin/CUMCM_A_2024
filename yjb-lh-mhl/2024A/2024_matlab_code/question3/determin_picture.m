function s = determin_picture(x_1, y_1, x1, y1)
    P = [x_1 y_1]; % 要检测的点即就是论文中的p点
    %%%%% 发生碰撞的板的4个坐标
    %%%%% 检测面积
    s1_x = [x_1 x1(1) x1(2)];
    s1_y = [y_1 y1(1) y1(2)];
    s1 = polyarea(s1_x, s1_y);

    s2_x = [x_1 x1(1) x1(3)];
    s2_y = [y_1 y1(1) y1(3)];
    s2 = polyarea(s2_x, s2_y);

    s3_x = [x_1 x1(3) x1(4)];
    s3_y = [y_1 y1(3) y1(4)];
    s3 = polyarea(s3_x, s3_y);

    s4_x = [x_1 x1(2) x1(4)];
    s4_y = [y_1 y1(2) y1(4)];
    s4 = polyarea(s4_x, s4_y);

    s5 = 0.66; % m²
    s = s4 + s3 + s2 + s1 -s5;
end
