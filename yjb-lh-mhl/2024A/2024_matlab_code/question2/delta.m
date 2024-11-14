function theta = delta(L, s)
    %用于求解两端点径向方向间初始夹角
    %输入:L 杆长
    %s 螺距宽
    %输出:两端点径向方向间初始夹角
    r = s * 16;
    syms theta
    f = @(theta) cos(theta) - ((s * (32 * pi - theta) / (2 * pi)) ^ 2 + r ^ 2 - L ^ 2) / (2 * (s * (32 * pi - theta) / (2 * pi)) * r);
    theta = fzero(f, 0.5);
end
