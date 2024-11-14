% 问题4结果分析
%% 初始化
clear
clc

%% 计算恒定的质心与转动惯量，存至全局变量之中
m1=4866;
R0=1;
h1=3;
h2=0.8;
global I1 rz1;
sgm=m1/(2*pi*R0*R0+2*pi*R0*h1+pi*R0*sqrt(R0^2+h2^2));
rz1=sgm*(2*pi*R0*h1+pi*R0*h2)/m1;
I11=(sgm*pi*R0^2)*(R0^2/2+(h1-rz1)^2+rz1^2);
I12=pi*R0^3*sgm*h1+2*pi*R0*sgm*((h1-rz1)^3/3+rz1^3/3);
I13=2*pi*sgm*R0* integral(@(z)((R0.*(h2-z)./h2).^2+(z+rz1).^2).*(h2-z)./h2,0,h2);
I1=I11+I12+I13;

%% main
[p1,p2]=aim([0,50000])

%% 目标函数
function [p1,p2]=aim(x)
    omg=1.9806;
    T=100*2*pi/omg;
    global C
    C(1)=x(1);
    C(2)=x(2);
    [~,n]=ode45('Q4_func1',0:0.01:100+T,[0,0,0,0,0,0,0,0]);
    p1=sum((C(1)*(abs(n(10000:end,4)-n(10000:end,3))).^2).*0.01)/T;
    p2=sum((C(2)*(abs(n(10000:end,8)-n(10000:end,7))).^2).*0.01)/T;
end