function dy = Q4_func1( t,y )
%% 参数赋值
m1=4866;
mf=1091.099;
R0=1;
h1=3;
h2=0.8;
m2=2433;
R1=0.5;
h3=0.5;
rou=1025;
g=9.8;
k=80000;
K=250000;
d0=0.5;
f=1760;
omg=1.9806;
lmd=528.5018;
mu=1655.909;
G=8890.7;
L=2140;
d=d0+(y(2)-y(1));
If=7142.493;

%% 微分方程组系数赋值
global C;
beta=C(1);
gamma=C(2);

%% 计算质心
global rz1
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz2=d+h3/2;
rz=(m1*rz1+m2*rz2)/(m1+m2);

%% 计算转动惯量
global I1
I2=m2*R0^2/4+m2*h3^2/12+m2*(d+h3/2)^2;

%% 计算相关参数
% delta=rz1/abs(rz1-rz);
delta=1;

%% 微分方程组描述
dy=zeros(8,1);
dy(1)=y(3);
dy(2)=y(4);
dy(3)=(f*cos(omg*t)-rou*g*y(1)*pi*R0*R0+k*(y(2)-y(1))+beta*(y(4)-y(3))-lmd*y(3))/(m1+mf);
dy(4)=(-m2*g-k*(y(2)-y(1))-beta*(y(4)-y(3)))/m2;
dy(5)=y(7);
dy(6)=y(8);
dy(7)=(delta*y(5)*G+delta*y(7)*mu+delta*L*cos(omg*t)+(y(7)-y(8))*gamma+(y(5)-y(6)*K))/(I1+If);
dy(8)=-((y(8)-y(7))*gamma+(y(6)-y(5))*K)/I2;
end