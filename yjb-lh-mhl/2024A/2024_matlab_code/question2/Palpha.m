function [alpha1, alpha2] = Palpha(r) % % % % % % 求两个 alpha 角度 / rad
    %%初始化参数
    L1 = 2.86; %龙头长度 单位：m
    L2 = 1.65; %龙身龙尾长度
    alpha1 = zeros(1, 223);
    alpha2 = alpha1;
    calpha1 = alpha1;
    calpha2 = alpha1;

    calpha1(1) = (L1 ^ 2 + r(1) .^ 2 - r(2) .^ 2) ./ (2 * L1 .* r(1));
    alpha1(1) = acos(calpha1(1));
    calpha2(1) = ((r(2)) .^ 2 + L1 ^ 2 - r(1) .^ 2) ./ (2 * r(2) .* L1);
    alpha2(1) = acos(calpha2(1));

    for i = 2:223
        calpha1(i) = (L2 ^ 2 + r(i) .^ 2 - r(i + 1) .^ 2) ./ (2 * L2 .* r(i));
        alpha1(i) = acos(calpha1(i));
        calpha2(i) = (r(i + 1) .^ 2 + L2 ^ 2 - r(i) .^ 2) ./ (2 * r(i + 1) .* L2);
        alpha2(i) = acos(calpha2(i));
    end

end
