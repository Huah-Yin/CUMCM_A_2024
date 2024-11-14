function xdot = que1_1(t, x, params) %parmas为要传递的参数
    mf = params.mf; %浮子质量:kg
    rf = params.rf; %浮子底半径：m
    mz = params.mz; %振子质量 (kg)
    rou = params.rou; %海水的密度 (kg/m3)
    g = params.g;
    k4 = params.k4; %垂荡兴波阻尼系数 (N·s/m)
    m0 = params.m0; %垂荡附加质量
    k2 = params.k2; %弹簧劲度系数
    A = params.A; %垂荡激励力振幅 (N)
    omega = params.omega; %入射波浪频率 (s-1)
    %%%%%%%%%question 1.2  detK为阻尼器变系数

    % x(1)=x1,x(2)=x1' ,x(3)=x2, x(4)=x2'
    % 浮子               振子
    detK = params.detK;

    xdot = [x(2); %浮子速度
            - (rou * g * pi * (rf ^ 2) * x(1) + k2 * (x(1) - x(3)) + detK * (x(2) - x(4)) + k4 * x(2) + A * cos(omega * t)) / (m0 + mf);
            x(4); %振子速度
            - (k2 * (x(3) - x(1)) + detK * (x(4) - x(2))) / mz
            ];
end
