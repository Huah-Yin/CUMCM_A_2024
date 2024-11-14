% 问题1参数估计
%% 初始化
clear
clc
close all
%% 设置炉温
global U_info
U_info=[ones(1,5)*175,195,235,ones(1,2)*255,ones(1,2)*25];
% 炉前,1-4,5,6,7,8,9,10-11,炉后
D_info=    [0.00014,ones(1,4)*0.00020,  0.00020,0.00020,0.00027,0.00020,0.00020,0.00015,0.00015];
sigma_info=[0.012,  ones(1,4)*0.021,    0.024,  0.030,  0.034,  0.025,0.025,0.010,0.010];
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
    for j=2:j_max-1
        data(k,j)=DD*(data(k-1,j+1)-2*data(k-1,j)+data(k-1,j-1))+data(k-1,j);
    end
    data(k,1)=-k_step*sigma*(data(k,1)-get_U(t_now*vv))+data(k-1,1);
    data(k,j_max)=-k_step*sigma*(data(k,j_max)-get_U(t_now*vv))+data(k-1,j_max);
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