function [r0, theta0, w0, v0] = new_prepare(t, dx, dt, n,v0)
    %%初始化参数
    L1 = 2.86; %龙头长度 单位：m
    L2 = 1.65; %龙身龙尾长度
    s = dx; %初始螺线宽距
    b = s / (2 * pi);
    delta_theta_1 = delta(L1, s, n);
    delta_theta_2 = delta(L2, s, n);
    t0 = 0:dt:t;
    t0 = t0';
    len_t = length(t0);
    v = zeros(len_t, 224); %定义线速度,一共224个把手
    theta = v;
    Omega = v; %定义角速度
    r = v; %定义极径
    v(:, 1) = v0;
    %%计算第一个点的角速度
    w1_funs = @(t, theta) v(1,1) ./ (b * sqrt(1 + (2 * n * pi - theta) .^ 2));
    theta(:, 1) = 2 * n * pi - dis_euler(w1_funs, 1, t / dt, 0, dt);
    r(:, 1) = b .* theta(:, 1);
    Omega(:, 1) = v(1,1) ./ sqrt(r(:, 1) .^ 2 + b ^ 2);

    %%计算第二个的角速度
    i = 2 * ones(1, 222);

    while (abs(theta(i(1), 1) - theta(1, 1)) < delta_theta_1)
        i (1) = i(1) + 1;

        if i(1) > len_t
            break;
        end

    end

    %刚性杆的约束条件
    calpha1 = @(t, theta_2) (L1 ^ 2 + r(t, 1) .^ 2 - (b * (2 * n * pi - theta_2)) .^ 2) ./ (2 * L1 .* r(t, 1));
    salpha1 = @(t, theta_2)sqrt(1 - calpha1(t, theta_2) .^ 2);
    calpha2 = @(t, theta_2)((b * (2 * n * pi - theta_2)) .^ 2 + L1 ^ 2 - r(t, 1) .^ 2) ./ (2 * b .* (2 * n * pi - theta_2) .* L1);
    salpha2 = @(t, theta_2)sqrt(1 - calpha2(t, theta_2) .^ 2);
    % 沿杆的分速度相等
    funs = @(t, theta_2) Omega(t, 1) .* (calpha1(t, theta_2) + theta(t, 1) .* salpha1(t, theta_2)) ./ (calpha2(t, theta_2) + (2 * n * pi - theta_2) .* salpha2(t, theta_2));
    y0 = 0;
    h = 0.1;
    theta(i(1):end, 2) = 2 * n * pi - dis_euler(funs, i(1), t / dt, y0, dt);
    r(i(1):end, 2) = b .* theta(i(1):end, 2);

    for k = i(1):length(t0)
        Omega(k, 2) = funs(k, 2 * n * pi - theta(k, 2));
    end

    v(i(1):end, 2) = Omega(i(1):end, 2) .* sqrt(r(i(1):end, 2) .^ 2 + b ^ 2);
    %计算第二个点以后的角速度
    flag = 0;

    for j = 3:224
        i(j - 1) = i(j - 2);

        while (abs(theta(i(j - 1) + 1, j - 1) - theta(i(j - 2) + 1, j - 1)) < delta_theta_2)
            i (j - 1) = i(j - 1) + 1;

            if i(j - 1) > (t / dt)
                flag = 1;
                break;
            end

        end

        if flag == 1
            break;
        end

        calpha1 = @(t, theta_2) (L2 ^ 2 + r(t, j - 1) .^ 2 - (b * (2 * n * pi - theta_2)) .^ 2) ./ (2 * L2 .* r(t, j - 1));
        salpha1 = @(t, theta_2)sqrt(1 - calpha1(t, theta_2) .^ 2);
        calpha2 = @(t, theta_2)((b * (2 * n * pi - theta_2)) .^ 2 + L2 ^ 2 - r(t, j - 1) .^ 2) ./ (2 * b .* (2 * n * pi - theta_2) .* L2);
        salpha2 = @(t, theta_2)sqrt(1 - calpha2(t, theta_2) .^ 2);
        funs = @(t, theta_2) Omega(t, j - 1) .* (calpha1(t, theta_2) + theta(t, j - 1) .* salpha1(t, theta_2)) ./ (calpha2(t, theta_2) + (2 * n * pi - theta_2) .* salpha2(t, theta_2));
        y0 = 0;
        h = 0.1;
        theta(i(j - 1):end, j) = 2 * n * pi - dis_euler(funs, i(j - 1), t / dt, y0, dt);
        r(i(j - 1):end, j) = b .* theta(i(j - 1):end, j);

        for k = i(j - 1):length(t0)
            Omega(k, j) = funs(k, 2 * n * pi - theta(k, j));
        end

        v(i(j - 1):end, j) = Omega(i(j - 1):end, j) .* sqrt(r(i(j - 1):end, j) .^ 2 + b ^ 2);
    end

  
    r0 = r;
    theta0 = theta;
    w0 = Omega;
    v0 = v;
end
