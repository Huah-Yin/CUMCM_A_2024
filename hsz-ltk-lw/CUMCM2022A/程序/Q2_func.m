function y = Q2_func(x)
%问题2遗传算法目标函数
y=0;
global C;
C(1)=x(1);
C(2)=x(2);
omg=1.7152;
T=50*2*pi/omg;
[t,n]=ode45('Q2_func1',[0:0.01:100+T],[0,0,0,0]);
u=n(:,4)-n(:,3);
su=0;
for i=10000:length(u)
    beta=C(1)*(abs(n(i,4)-n(i,3)))^C(2);
    su=su+beta*u(i)*u(i)*0.01;
end
su=su/T;
y=su;
end