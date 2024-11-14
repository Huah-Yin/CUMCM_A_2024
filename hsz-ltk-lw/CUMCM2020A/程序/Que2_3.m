% 问题2结果分析
%% 初始化
clear
clc
close all
%% 设置炉温
global U_info
U_info=[ones(1,5)*182,203,237,ones(1,2)*254,ones(1,2)*25];
% 炉前,1-4,5,6,7,8,9,10-11,炉后
D_info=    [0.00014,ones(1,4)*0.00020,  0.00020,0.00020,0.00027,0.00020,0.00020,0.00015,0.00015];
sigma_info=[0.012,  ones(1,4)*0.021,    0.024,  0.030,  0.034,  0.025,0.025,0.010,0.010];
len=30.5;
gap=5;
global x0 y0 D_ls sigma_ls
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

plt_init()
for v=81.5691 %65:10:100
    res(v)=constraint(v);
end
legend(["约束2(1)","约束2(2)","约束3","约束4(1)","约束4(2)","过炉速度"+[65:10:100]+"cm/s"])
xlabel("时间/s")
ylabel("温度/℃")

function res=constraint(v)
% 约束条件
% 输入速度v(cm/s);输出1-符合约束,0-不符合约束
    global D_ls sigma_ls
    d=0.15;     % 厚度(mm)
    j_step=0.005;   % 厚度步长(mm)
    k_step=0.005;     % 时间步长(s)
    vv=v/60;    % 传送带前进速度(cm/s)[不要修改参数]
    data=zeros(floor(435/(v/60)/k_step)+1,floor(d/j_step)+1);     % (k,j)->(时间s,空间cm);首行/首列为零时刻/位置
    [k_max,j_max]=size(data);
    data(1,:)=25;   % 初始条件
    key=ones(k_max,1)*25; % 中心温度记录
    count1=0;   % 约束2计时器
    count2=0;   % 约束3计时器
    for k=2:k_max
        t_now=k*k_step; % 当前时间(s)
        DD=k_step*get_D(t_now*vv,D_ls)/(j_step*j_step);
        sigma=get_sigma(t_now*vv,sigma_ls);
        for j=2:j_max-1
            data(k,j)=DD*(data(k-1,j+1)-2*data(k-1,j)+data(k-1,j-1))+data(k-1,j);
        end
        data(k,1)=-k_step*sigma*(data(k-1,1)-get_U(t_now*vv))+data(k-1,1);
        data(k,j_max)=-k_step*sigma*(data(k-1,j_max)-get_U(t_now*vv))+data(k-1,j_max);
        key(k)=data(k,d/2/j_step);
        
        if key(k)-key(k-1)>3*k_step || key(k)-key(k-1)<-3*k_step    % 约束(1)
            res=0;
            fprintf("v=%f,约束(1)跳出\n",v)
            return
        end
        if key(k)-key(k-1)>0 && key(k)>150 && key(k)<190
            count1=count1+1;
        elseif key(k)>217
            count2=count2+1;
        end
        if count1*k_step>120 || count2*k_step>90 || key(k)>250   % 约束(2)(3)(4)上界
%             plt(k,key,t_now,k_step)
            plot(0:k_step:t_now-k_step,key(1:k),"k")
            fprintf("v=%f,约束(234)上界跳出\n",v)
            res=0;
            return
        end
    end
    if count1*k_step<60 || count2*k_step<40 || max(key)<240 % 约束(2)(3)(4)下界
        plot(0:k_step:t_now-k_step,key(1:k))
        fprintf("v=%f,约束(234)下界跳出\n",v)
        res=0;
        return
    end
    fprintf("v=%f,符合条件\n",v)
    plot(0:k_step:t_now-k_step,key(1:k),"k")
    res=1;
end

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

function plt_init()
% 绘图
    figure
    hold on
    plot([0,450],[150,150],"r--", ...
        [0,450],[190,190],"r--", ...
        [0,450],[217,217],"g--", ...
        [0,450],[240,240],"b--", ...
        [0,450],[250,250],"b--")
end