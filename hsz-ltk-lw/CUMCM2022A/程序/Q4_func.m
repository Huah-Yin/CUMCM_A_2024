function y = Q4_func(x)
%问题4遗传算法目标函数
y=0;
global C;
C(1)=x(1);
C(2)=x(2);
omg=1.9806;
T=50*2*pi/omg;
[t,n]=ode45('Q4_func1',0:0.01:100+T,[0,0,0,0,0,0,0,0]);
u=n(:,4)-n(:,3);
w=n(:,8)-n(:,7);
su=0;
beta=C(1);
gamma=C(2);
for i=10000:length(u)
    su=su+(beta*u(i)*u(i)+gamma*w(i)*w(i))*0.01;
end
su=su/T;
y=su;
end