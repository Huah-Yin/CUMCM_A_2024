% 问题3与问题2对照
%% 初始化
clear
clc
close all

%% 设置炉温
U_info2=[ones(1,5)*173,198,230,ones(1,2)*257,ones(1,2)*25];
U_info3=[ones(1,5)*176.2229,190.3380,230.9501,ones(1,2)*264.6106,ones(1,2)*25];
len=30.5;
gap=5;
x0=[0,25];
y2=[25];
y3=[25];
for i=1:11
    x0(end+1)=x0(end)+len;
    x0(end+1)=x0(end)+gap;
    y2(end+1)=U_info2(i);
    y2(end+1)=U_info2(i);
    y3(end+1)=U_info3(i);
    y3(end+1)=U_info3(i);
end
x0(end)=435.5;
y2(end+1)=25;
y3(end+1)=25;

%% 绘制炉内温度图
x_fig=0:0.5:435.5;
y2_fig=[];
y3_fig=[];
for x=x_fig
    y2_fig(end+1)=interp1(x0,y2,x);
    y3_fig(end+1)=interp1(x0,y3,x);
end
plot(x_fig,y2_fig,"k",x_fig,y3_fig,"r")
xlabel("炉内位置/cm")
ylabel("温度/℃")
xlim([0,435.5])
legend("问题2","问题3")

% 问题3炉温曲线绘图
%% 初始化
figure
hold on
U_ls=[U_info2;U_info3];
for idx=[1,2]
    %% 设置炉温
    global U_info
    U_info=U_ls(idx,:);
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

    d=0.15;     % 厚度(mm)
    j_step=0.005;   % 厚度步长(mm)
    k_step=0.005;     % 时间步长(s)
    v=[81.5691,90.7978];   % 传送带前进速度(cm/min)
    vv=v(idx)/60;    % 传送带前进速度(cm/s)[不要修改参数]
    data=zeros(floor(435/vv/k_step)+1,floor(d/j_step)+1);     % (k,j)->(时间s,空间cm);首行/首列为零时刻/位置
    [k_max,j_max]=size(data);
    data(1,:)=25;   % 初始条件
    for k=2:k_max
        t_now=k*k_step; % 当前时间(s)
        DD=k_step*get_D(t_now*vv,D_ls)/(j_step*j_step);
        sigma=get_sigma(t_now*vv,sigma_ls);
        for j=2:j_max-1
            data(k,j)=DD*(data(k-1,j+1)-2*data(k-1,j)+data(k-1,j-1))+data(k-1,j);
        end
        data(k,1)=-k_step*sigma*(data(k-1,1)-get_U(t_now*vv))+data(k-1,1);
        data(k,j_max)=-k_step*sigma*(data(k-1,j_max)-get_U(t_now*vv))+data(k-1,j_max);
    end
    plot(0:k_step:435/vv,data(:,d/2/j_step))
    xlabel("时间/s")
    ylabel("温度/℃")
end
plot([0,350],[217,217],'r--')
legend("问题2","问题3","217℃线")


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