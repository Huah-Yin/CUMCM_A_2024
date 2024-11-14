%% 初始化
clear
close all
clc

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
d0=0.5;
f=3640;
omg=1.7152;
lmd=683.4558;
mu=654.3383;
beta=10000;
gamma=1000;
syms d;

%% 计算质心
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz1=sgm*(pi*R0*h1^2+2*pi*R0/h2*(1+R0/h2)*(-(h2^3)/6)+h1*pi*R0^2)/m1;
rz2=d+h3/2;
rz=(m1*rz1+m2*rz2)/(m1+m2);

%% 计算转动惯量
I11=(sgm*pi*R0^2)*(R0^2/2+(h1-rz1)^2+rz1^2);
I12=pi*R0^3*sgm*h1+2*pi*R0*sgm*((h1-rz1)^3/3+rz1^3/3);
I13=2*pi*sgm*R0* integral(@(z)((R0.*(h2-z)./h2).^2+(z+rz1).^2).*(h2-z)./h2,0,h2);
I1=I11+I12+I13;
I2=m2*R0^2/4+m2*h3^2/12+m2*(d+h3/2)^2;
t0=0;
tf=40*2*pi/omg;

%% 求解微分方程组、作图、保存数据
[t,n]=ode45('Q3_func',t0:0.2:tf,[0,0,0,0,0,0,0,0]);
figure
subplot(2,2,1)
hold on
plot(t(:),n(:,1),'r-');
plot(t(:),n(:,2),'b-');
legend('浮子位移','振子位移')
xlabel('时间(s)');
ylabel('位移(m)');
title('垂荡位移')
hold off
subplot(2,2,2)
hold on
plot(t(:),n(:,3),'r-');
plot(t(:),n(:,4),'b-');
legend('浮子速度','振子速度')
xlabel('时间(s)');
ylabel('速度(m/s)');
title('垂荡速度')
hold off
subplot(2,2,3)
hold on
plot(t(:),n(:,5),'r-');
plot(t(:),n(:,6),'b-');
legend('浮子角位移','振子角位移')
xlabel('时间(s)');
ylabel('角位移(rad)');
title('纵摇角位移')
hold off
subplot(2,2,4)
hold on
plot(t(:),n(:,7),'r-');
plot(t(:),n(:,8),'b-');
legend('浮子角速度','振子角速度')
xlabel('时间(s)');
ylabel('角速度(rad/s)');
title('纵摇角速度')
hold off
A=t;
A(:,2)=n(:,1);
A(:,3)=n(:,3);
A(:,4)=n(:,5);
A(:,5)=n(:,7);
A(:,6)=n(:,2);
A(:,7)=n(:,4);
A(:,8)=n(:,6);
A(:,9)=n(:,8);
xlswrite('问题3结论.xlsx',A,1);
ind=1;
for i=1:length(A(:,1))
    if ismember(A(i,1),[10 20 40 60 100])
        C(ind,1:9)=A(i,1:9);
        ind=ind+1;
    end
end