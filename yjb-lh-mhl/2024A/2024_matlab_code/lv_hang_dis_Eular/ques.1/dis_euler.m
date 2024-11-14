function y=dis_euler(f, a, b, y0)
h =0.1;
s = b - a+1;
Y = zeros(1,s+1);
Y(1) = y0;
for k = 1:s
    Y(k+1) = Y(k) + h * f(a+k-1, Y(k));
end
y = Y';%转置
end
