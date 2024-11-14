% 获得P与rou的离散关系
%% 初始化
clear
clc
close all

%% 读附件3
data3=xlsread("附件3-弹性模量与压力.xlsx");
x=data3(:,1);
y=data3(:,2);

%% 左侧循环积分（前缀和）
h=0.5;  % 步长
down = 100;
c=0;
result = [];
for i=1:length(y)-1
    c=(1/y(i)+1/y(i+1))/2*(x(2)-x(1))+c ; 
    result(end+1) = c;
end
%常数追加
det=result(find(x==down));
result=result-det;

%% 与右侧积分相等，反解rou
rou=exp(result+log(0.85));

%% 绘图
plot(x(1:end-1),rou,"o-")
xlabel("压强P")
ylabel("密度ρ")
title("密度与压强的关系图")

%% 导出离散数据
dataout={"压强(MPa)","密度(mg/mm^3)"};
for ind=1:length(result)
    dataout{ind+1,1}=x(ind);
    dataout{ind+1,2}=rou(ind);
end
xlswrite("压强与密度的关系.xlsx",dataout,1);

%% 拟合
p=polyfit(x(1:end-1),rou,1);
func=@(x)p(1).*x+p(2);

%% 检验
TSS=sum((rou-mean(rou)).^2);
RSS=sum((func(x(1:end-1)')-rou).^2);
R2=1-RSS/TSS

