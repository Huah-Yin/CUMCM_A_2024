function [alpha1, alpha2] = alpha(r, b, theta, L1, L2, i)
    calpha1 = (L1 ^ 2 + r(i(1):end, 1) .^ 2 - r(i(1):end, 2) .^ 2) ./ (2 * L1 .* r(i(1):end, 1));
    alpha1(:, 1) = acos(calpha1);
    calpha2 = (r(i(1):end, 2) .^ 2 + L1 ^ 2 - r(i(1):end, 1) .^ 2) ./ (2 .* r(i(1):end) .* L1);
    alpha2(:, 1) = acos(calpha2);

    for j = 2:223
        calpha3 = (L2 ^ 2 + r(i(n):end, n) .^ 2 - (b * (32 * pi - theta(i(n + 1:end, n + 1)) .^ 2))) ./ (2 * L2 .* r(i(n):end, 1));
        alpha1(:, j) = acos(calpha3);
        calpha4 = ((b * (32 * pi - theta(:, j + 1))) .^ 2 + L2 ^ 2 - r(:, j) .^ 2) ./ (2 * b .* (32 * pi - theta(:, j + 1)) .* L2);
        alpha2(:, j) = acos(calpha4);
    end

end
