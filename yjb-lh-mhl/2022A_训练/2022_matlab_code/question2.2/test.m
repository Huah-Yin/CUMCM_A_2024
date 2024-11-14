function Y = test(X)
    d = length(X); % dimensions of the problem

    % Parameters of the problem
    a = 20;
    b = 0.2;
    c = 2 * pi;

    A1 = 0;
    A2 = 0;

    for i = 1:d
        xi = X(i); % variable xi
        A1 = A1 + xi ^ 2;
        A2 = A2 + cos(c * xi);
    end

    Y = -a * exp(-b * sqrt(A1 / d)) - exp(A2 / d) + a + exp(1);
end
