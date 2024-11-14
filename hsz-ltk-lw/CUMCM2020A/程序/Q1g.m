% 问题1参数估计
%% 初始化
clear
clc
close all
%% 设置炉温
global U_info
U_info=[ones(1,5)*175,195,235,ones(1,2)*255,ones(1,2)*25];
% 炉前,1,2-4,5,6,7,8,9,10-11,炉后
D_info=    1e-5*[0.1,0.24,ones(1,3)*0.27,  0.20,0.23,0.30,0.30,0.32,0.12,0.3];
sigma_info=[0.2, 2.2 ,ones(1,3)*3.3,    2.4,  3.0,  3.4,  2.5,2.5,0.4,1.2];
len=30.5;
gap=5;
global x0 y0
x0=[0,25];
y0=[25];
D_ls=[D_info(1)];
sigma_ls=[sigma_info(1)];
for i=1:11
    x0(end+1)=x0(end)+len;
    x0(end+1)=x0(end)+gap;
    y0(end+1)=U_info(i);
    y0(end+1)=U_info(i);
    D_ls(end+1)=D_info(i+1);
    D_ls(end+1)=D_info(i+1);
    sigma_ls(end+1)=sigma_info(i+1);
    sigma_ls(end+1)=sigma_info(i+1);
end
x0(end)=435.5;
y0(end+1)=25;
D_ls(end+1)=D_ls(end);
sigma_ls(end+1)=sigma_info(end);

%% 绘制炉内温度图
% x_fig=0:0.2:435.5;
% y_fig=[];
% D_fig=[];
% sigma_fig=[];
% for x=x_fig
%     y_fig(end+1)=get_U(x);
%     D_fig(end+1)=get_D(x,D_ls);
%     sigma_fig(end+1)=get_sigma(x,sigma_ls);
% end
% plot(x_fig,y_fig,"k")
% xlabel("炉内位置/cm")
% ylabel("温度/℃")
% xlim([0,435.5])
% figure
% plot(x_fig,D_fig,"b")
% figure
% plot(x_fig,sigma_fig,"g")

%% 导入数据
aim=xlsread("附件.xlsx");

%% 求解炉温曲线
d=0.15;     % 厚度(mm)
j_step=0.005;   % 厚度步长(mm)
k_step=0.005;     % 时间步长(s)
v=70;   % 传送带前进速度(cm/min)
vv=v/60;    % 传送带前进速度(cm/s)[不要修改参数]
data=zeros(floor(435/(v/60)/k_step)+1,floor(d/j_step)+1);     % (k,j)->(时间s,空间cm);首行/首列为零时刻/位置
[k_max,j_max]=size(data);
data(1,:)=25;   % 初始条件
for k=2:k_max
    t_now=k*k_step; % 当前时间(s)
    DD=k_step*get_D(t_now*vv,D_ls)/(j_step*j_step);
    sigma=get_sigma(t_now*vv,sigma_ls);
    r=DD*k_step/(j_step*j_step);
    A=eye(j_max)*2*r+1;
    A(1,1)=1+j_step*sigma;
    A(end,end)=1+j_step*sigma;
    A(1,2)=-1;
    A(end,end-1)=-1;
    for iidx=2:j_max-1
        A(iidx,iidx-1)=-r;
        A(iidx,iidx+1)=-r;
    end
    B=data(k-1,:)';
    B(1)=j_step*sigma*get_U(t_now*vv);
    B(end)=j_step*sigma*get_U(t_now*vv);
    data(k,:)=Chase_method(A,B);
end
%% 计算误差平方和
err=sum(abs(aim(1:end-1,2)-data(aim(1,1)/k_step:0.5/k_step:min(aim(end,1)/k_step,k_max),d/2/j_step))./aim(1:end-1,2));
errmean=err/(length(aim)-1);
[X,Y]=meshgrid(0:k_step:435/(v/60),0:j_step:d);
mesh(X,Y,data')
xlabel("时间/s")
ylabel("深度/mm")
zlabel("温度/℃")
colorbar()
figure
plot(0:k_step:435/(v/60),data(:,d/2/j_step),"b",aim(:,1),aim(:,2),"g",0:k_step:435/(v/60),data(:,1),"k--",0:k_step:435/(v/60),get_U([0:k_step:435/(v/60)].*vv),"r--")
legend("当前中心炉温曲线","标准炉温曲线","电路板表面温度","炉内空气温度")
xlabel("时间/s")
ylabel("温度/℃")
% figure
% plot(0:j_step:d,data(400,:))

function U=get_U(x)
% 输入x-炉内位置(cm)，输出对应炉内温度
    global x0 y0
    U=interp1(x0,y0,x);
end

function D=get_D(x,D_ls)
% 输入x-炉内位置(cm)，输出对应炉内温度
    global x0
    D=interp1(x0,D_ls,x);
end

function sigma=get_sigma(x,sigma_ls)
% 输入x-炉内位置(cm)，输出对应炉内温度
    global x0
    sigma=interp1(x0,sigma_ls,x);
end

function [ x ] = Chase_method( A, b )
%Chase method 追赶法求三对角矩阵的解
%   A为三对角矩阵的系数，b为等式右端的常数项，返回值x即为最终的解
%   注：A尽量为方阵，b一定要为列向量
%% 求追赶法所需L及U
T = A;
for i = 2 : size(T,1)
    T(i,i-1) = T(i,i-1)/T(i-1,i-1);
    T(i,i) = T(i,i) - T(i-1,i) * T(i,i-1);
end
L = zeros(size(T));
L(logical(eye(size(T)))) = 1;   %对角线赋值1
for i = 2:size(T,1)
    for j = i-1:size(T,1)
        L(i,j) = T(i,j);
        break;
    end
end
U = zeros(size(T));
U(logical(eye(size(T)))) = T(logical(eye(size(T))));
for i = 1:size(T,1)
    for j = i+1:size(T,1)
        U(i,j) = T(i,j);
        break;
    end
end
%% 利用matlab解矩阵方程的遍历直接求解
y = L\b;
x = U\y;
end