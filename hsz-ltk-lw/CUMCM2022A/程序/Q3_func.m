function dy = Q3_func( t,y )
% 问题3微分方程组

%% 参数赋值
m1=4866;
mf=1028.876;
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
f=3640;
omg=1.7152;
lmd=683.4558;
mu=654.3383;
beta=10000;
gamma=1000;
G=8890.7;
L=1690;
d=d0+(y(2)-y(1));
If=7001.914;

%% 计算质心
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz1=sgm*(2*pi*R0*h1+pi*R0*h2)/m1;
rz2=d+h3/2;
rz=(m1*rz1+m2*rz2)/(m1+m2);

%% 计算转动惯量
I11=(sgm*pi*R0^2)*(R0^2/2+(h1-rz1)^2+rz1^2);
I12=pi*R0^3*sgm*h1+2*pi*R0*sgm*((h1-rz1)^3/3+rz1^3/3);
I13=2*pi*sgm*R0* integral(@(z)((R0.*(h2-z)./h2).^2+(z+rz1).^2).*(h2-z)./h2,0,h2);
I1=I11+I12+I13;
I2=m2*R0^2/4+m2*h3^2/12+m2*(d+h3/2)^2;

%% 计算相关参数
delta=abs(rz1-rz)/rz1;

%% 微分方程组描述
dy=zeros(8,1);
dy(1)=y(3);
dy(2)=y(4);
dy(3)=(f*cos(omg*t)-rou*g*y(1)*pi*R0*R0+k*(y(2)-y(1))+beta*(y(4)-y(3))-lmd*y(3))/(m1+mf);
dy(4)=(-k*(y(2)-y(1))-beta*(y(4)-y(3)))/m2;
dy(5)=y(7);
dy(6)=y(8);
dy(7)=(delta*y(5)*G+delta*y(7)*mu+delta*L*cos(omg*t)+(y(7)-y(8))*gamma+(y(5)-y(6)*K))/(I1+If);
dy(8)=-((y(8)-y(7))*gamma+(y(6)-y(5))*K)/I2;
end