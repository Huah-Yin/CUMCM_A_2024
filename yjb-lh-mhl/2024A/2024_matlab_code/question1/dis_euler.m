function y = dis_euler(f, a, b, y0) % % % % % % %  欧拉法
    dt = 0.1;
    s = b - a + 1;
    Y = zeros(1, s + 1);
    Y(1) = y0;

    for k = 1:s
        y_pred = Y(k) + dt * f(a + k - 1, Y(k));

        % 修正步骤（使用预测值进行修正）
        Y(k + 1) = Y(k) + (dt / 2) * (f(a + k - 1, Y(k)) + f(a + k, y_pred));
        % Y(k + 1) = Y(k) + h * f(a + k - 1, Y(k));
    end

    y = Y';
end
