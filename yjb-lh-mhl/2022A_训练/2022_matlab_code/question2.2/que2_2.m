function ydot = que2_2(t, y, params, x)
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
    %%%%函数
    % y(1)=x1,y(2)=x1' ,y(3)=x2, y(4)=x2'
    % 浮子               振子
    detK = x(1) * power(abs(y(2) - y(4)),x(2)); %x(1)为系数k,x(2)为系数:阿尔法
    ydot = [y(2);
            - (rou * g * pi * (rf ^ 2) * y(1) + k2 * (y(1) - y(3)) + detK * (y(2) - y(4)) + k4 * y(2) + A * cos(omega * t)) / (m0 + mf);
            y(4);
            - (k2 * (y(3) - y(1)) + detK * (y(4) - y(2))) / mz
            ];
end
