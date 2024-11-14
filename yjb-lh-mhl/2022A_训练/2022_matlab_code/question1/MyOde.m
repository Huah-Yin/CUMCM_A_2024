function [x, t] = MyOde(f, x0, t) %采取rk4
    nt = numel(t);
    nx = numel(x0);

    x = nan(nx, nt);

    x(:, 1) = x0;

    for k = 1:nt - 1
        k1 = dt * f(t(k), x(:, k));
        k2 = dt * f(t(k) + dt / 2, x(:, k) + k1 / 2);
        k3 = dt * f(t(k) + dt / 2, x(:, k) + k2 / 2);
        k4 = dt * f(t(k) + dt, x(:, k) + k3);

        dx = (k1 + k2 * 2 + 2 * k3 + k4) / 6;
        x(:, k + 1) = x(:, k) + dx;
    end

end
