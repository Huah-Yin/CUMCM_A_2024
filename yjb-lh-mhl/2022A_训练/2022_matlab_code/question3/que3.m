function xdot = que3(t, x, params) %parmas为要传递的参数
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
    %%%纵摇参数
    L1 = params.L1;
    L2 = params.L2;
    L3 = params.L3;
    L4 = params.L4;
    L = params.L;
    Jc = params.Jc;
    o = params.o;
    x0 = mz * g / k2;
    %%%%%%%%%question 1.2  detK为阻尼器变系数

    % x(1)=x1,x(2)=x1' ,x(3)=x2, x(4)=x2'
    % 浮子               振子
    % x(5)=theta1 x(6)=theta1' x(7)=theta2  x(8) =theta2'
    %角加速度1                          角加速度2

    detK = params.detK;

    xdot = [x(2); %浮子速度
            - (rou * g * pi * (rf ^ 2) * x(1) + k2 * (x(1) - x(3)) + detK * (x(2) - x(4)) + k4 * x(2) + A * cos(omega * t)) / (m0 + mf);
            x(4); %振子速度
            - (k2 * (x(3) - x(1)) + detK * (x(4) - x(2))) / mz;

            x(6); %角速度1
            (-L1 * x(5) * o / (1/3 * (1/4 + x(3) - x(1) + x0 + 1.236)) - L2 * (x(5) - x(7)) - L3 * (x(6) - x(8)) - L4 * o / (1/3 * (1/4 + x(3) - x(1) + x0 + 1.236)) * x(6) + o / (1/3 * (1/4 + x(3) - x(1) + x0 + 1.236)) * L * cos(omega * t)) / (Jc + 0.5 * mf * (1.236 + x(3) - x(1) + x0 + 0.25) .^ 2/9);
            x(8); %角速度2
            (L2 * (x(5) - x(7)) + L3 * (x(6) - x(8))) / (mz * 0.5/3 * (3 * (x(3) - x(1) + x0) ^ 2 + 3 * 0.5 * (x(3) - x(1) + x0) + 0.5 * 0.5));
            ];
end
