function dy = Q2_func1( t,y )
%% 参数赋值
m1=4866;
mf=1165.992;
R0=1;
h1=3;
h2=0.8;
m2=2433;
R1=0.5;
h3=0.5;
rou=1025;
g=9.8;
k=80000;
d0=0.5;
f=4890;
omg=2.2143;
lmd=167.8395;

%% 微分方程组系数赋值
global C;
beta=C(1)*(abs(y(4)-y(3)))^C(2);

%% 微分方程组描述
dy=zeros(4,1);
dy(1)=y(3);
dy(2)=y(4);
dy(3)=(f*cos(omg*t)-rou*g*y(1)*pi*R0*R0+k*(y(2)-y(1))+beta*(y(4)-y(3))-lmd*y(3))/(m1+mf);
dy(4)=(-k*(y(2)-y(1))-beta*(y(4)-y(3)))/m2;
end