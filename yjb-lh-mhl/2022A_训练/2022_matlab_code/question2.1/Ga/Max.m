function P_max = Max(params_x)
    params.mf = 4866; %浮子质量:kg
    params.rf = 1; %浮子底半径：m
    params.mz = 2433; %振子质量 (kg)
    params.rou = 1025; %海水的密度 (kg/m3)
    params.g = 9.8;
    params.k4 = 167.8395; %垂荡兴波阻尼系数 (N·s/m)
    params.m0 = 1165.992; %垂荡附加质量
    params.k2 = 80000; %弹簧劲度系数
    params.A = 4890; %垂荡激励力振幅 (N)
    params.omega = 2.2143; %入射波浪频率 (s-1)
    %%%步长
    dt = 0.2;
    t = 0:dt:200;
    y0 = [0 0 0 0];
    %%%%函数
    % y(1)=x1,y(2)=x1' ,y(3)=x2, y(4)=x2'
    % 浮子               振子
    f = @(t, y) que2_2(t, y, params, params_x);
    %%求解
    [T, X] = ode45(f, t, y0);

    t1 = 100; %积分开始区间
    t_end = 180; %30个波浪周期

    idx_left = round(t1 / dt) + 1;
    idx_right = round(t_end / dt) + 1;
    % 检查并调整索引范围
    idx_left = max(idx_left, 1);
    idx_right = min(idx_right, length(T));

    tspan = idx_left:idx_right;

    v1 = X(tspan, 2);
    v2 = X(tspan, 4);

    T_span = T(tspan);

    diff_2 = params_x * power(abs(v1 - v2), 2);

    P_max = -((1 / (t_end - t1)) * trapz(T_span, diff_2)); %trapz 自带/2
end
