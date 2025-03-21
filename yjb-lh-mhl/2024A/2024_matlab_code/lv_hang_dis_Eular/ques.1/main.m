clc,clear,close all;
%%初始化参数
L1 = 2.86;%龙头长度 单位：m
L2 = 1.65;%龙身龙尾长度
s = 0.55;%初始螺线宽距
b = s/(2*pi);
delta_theta_1 = delta(L1,s);
delta_theta_2 = delta(L2,s);
t0 = 0:0.1:300;
t0 = t0';
v = zeros(length(t0),224);%定义线速度
 theta = v;
w =v;%定义角速度
r = v;%定义极径
v(:,1) = 1;
%%计算第一个点的角速度

w1_funs = @(t,theta) 1./(b*sqrt(1+(32*pi-theta).^2));
theta(:,1) = 32*pi-dis_euler(w1_funs,1,3000,0);
r(:,1) = b.*theta(:,1);
w(:,1) = 1./sqrt(r(:,1).^2+b^2);

%%计算第二个的角速度
i = 2*ones(1,222);
while(delta_theta_1-abs(theta(i(1),1)-theta(1,1))>0.005)
    i (1) = i(1)+1;
end
%刚性杆的约束条件
calpha1 =@(t,theta_2) (L1^2+r(t,1).^2-(b*(32*pi-theta_2)).^2)./(2*L1.*r(t,1));
salpha1 = @(t,theta_2)sqrt(1-calpha1(t,theta_2).^2);
calpha2 = @(t,theta_2)((b*(32*pi-theta_2)).^2+L1^2-r(t,1).^2)./(2*b.*(32*pi-theta_2).*L1);
salpha2 =@(t,theta_2)sqrt(1-calpha2(t,theta_2).^2);
% 沿杆方向速度相等
funs = @(t,theta_2) w(t,1).*(calpha1(t,theta_2)+theta(t,1).*salpha1(t,theta_2))./(calpha2(t,theta_2)+(32*pi-theta_2).*salpha2(t,theta_2));
y0 =0;
h = 0.1;
theta(i(1):end,2) = 32*pi-dis_euler(funs,i(1),3000,y0);
r(i(1):end,2) = b.*theta(i(1):end,2);
for k = i(1):length(t0)
    w(k,2) = funs(k,32*pi-theta(k,2));
end
v(i(1):end,2) = w(i(1):end,2).*sqrt(r(i(1):end,2).^2+b^2);
%计算第二个点以后的角速度
for j = 3:224
    i(j-1) = i(j-2);
    while(delta_theta_2-abs(theta(i(j-1)+1,j-1)-theta(i(j-2)+1,j-1))>0.005)
    i (j-1) = i(j-1)+1;
    end
calpha1 =@(t,theta_2) (L1^2+r(t,j-1).^2-(b*(32*pi-theta_2)).^2)./(2*L1.*r(t,j-1));
salpha1 = @(t,theta_2)sqrt(1-calpha1(t,theta_2).^2);
calpha2 = @(t,theta_2)((b*(32*pi-theta_2)).^2+L1^2-r(t,j-1).^2)./(2*b.*(32*pi-theta_2).*L1);
salpha2 =@(t,theta_2)sqrt(1-calpha2(t,theta_2).^2);
funs = @(t,theta_2) w(t,j-1).*(calpha1(t,theta_2)+theta(t,j-1).*salpha1(t,theta_2))./(calpha2(t,theta_2)+(32*pi-theta_2).*salpha2(t,theta_2));
y0 =0;
h = 0.1;
theta(i(j-1):end,j) = 32*pi-dis_euler(funs,i(j-1),3000,y0);
r(i(j-1):end,j) = b.*theta(i(j-1):end,j);
for k = i(j-1):length(t0)
    w(k,j) = funs(k,32*pi-theta(k,j));
end
v(i(j-1):end,j) = w(i(j-1):end,j).*sqrt(r(i(j-1):end,j).^2+b^2);
end
%%转化为直角坐标
[x,y] = pol2cart(theta,r);
y(1,1) = 0;
for n = 1:223
    y(i(n),n+1) = 0;
end