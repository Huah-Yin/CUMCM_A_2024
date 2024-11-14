%% 初始化
clear
clc
close all

%% 读文件
data=xlsread("附件1-凸轮边缘曲线.xlsx");

%% 拟合
f=fit(data(:,1),data(:,2),"fourier1");

%% 绘图
polarplot(data(:,1),data(:,2),"o",data(:,1),f(data(:,1)),"-")
legend("实际值","拟合曲线")

%% 检验
rou=data(:,2);
theta=data(:,1);
TSS=sum((rou-mean(rou)).^2);
RSS=sum((f(theta)-rou).^2);
R2=1-RSS/TSS