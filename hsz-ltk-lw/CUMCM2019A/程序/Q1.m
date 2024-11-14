%% 初始化
clear
clc
close all

%% 常量定义
global C A V0 get_rou get_P Qout Qin dt tend
C=0.85;
A=pi*(1.4/2)^2;
V0=pi*500*(10/2)^2;

%% rou与P的关系
get_rou=@(p)3.910284764226142e-04.*p+0.808708414066112;
get_P=@(rou)(rou-0.808708414066112)/3.910284764226142e-04;

%% 流出流量函数
Qout=@(t)(100*t).*(t<=0.2)+(20).*(t>0.2 && t<2.2)+(-100*t+240).*(t>=2.2 && t<=2.4)+0.*(t>2.4);

%% 流入流量函数
Qin=@(P2,P1)(C*A*(2*(P2-P1)/get_rou(160))^(0.5));

%% 循环优化Tin
dt=0.001;   % 差分步长，单位ms
tend=30000;    % 终止时间，单位ms
t1=0.01;     % 大步长
t2=0.0001;  % 小步长
%% 大步长搜索
% p=[];
% for Tin=10:t1:5000
%     p(end+1)=get_dPmax(Tin);
% end
% plot(10:t1:5000,p)

%% 小步长搜索
% p=[];
% for Tin=10.2:t2:10.4
%     p(end+1)=get_dPmax(Tin);
% end
% x=10.2:t2:10.4;
% plot(x,p,"b+-")
% xlabel("Tin(ms)")
% ylabel("Z(MPa)")
% [~,min_ind]=min(p);

% fprintf("最优周期:%f\n",x(min_ind))
dPmax=get_dPmax(10.2876);
%% 函数计算压强变化峰值
function dPmax = get_dPmax(Tin)
    global C A V0 get_rou get_P Qout Qin dt tend
    %Tin-流入周期
    Tout=100;
    P0=100; % 初始压强
    P2=160; % 进入压强
    m=get_rou(P0)*V0;  % 实时质量，赋初值为初始质量
    tin=1;  % 流入计时器
    tout=1; % 流出计时器
    ind=1;
    P=P0;   % 实时压强存储
    p_ls=P;
    dPmax=0;    % 压强差最大值
    for t=1:dt:tend
        tic
        %% 计时器重置
        if tin>Tin
            tin=tin-Tin;
        end
        if tout>Tout
            tout=tout-Tout;
        end
        %% 流入
        if tin<Tin-10
            vin=Qin(P2,P);
        else
            vin=0;
        end
        %% 流出
        vout=Qout(tout);
        %% 质量更新
        m=m+get_rou(160)*vin*dt-get_rou(P)*vout*dt;
        %% 压强变化计算
        rou_now=m/V0;
        P=get_P(rou_now);
        p_ls(end+1)=P;
        %% 计时器更新
        tin=tin+dt;
        tout=tout+dt;
        %% 越界判断
        if P<=0 || P>=200
            dPmax=inf;
            return
        end
        if dPmax<abs(P-P0)
            dPmax=abs(P-P0);
        end
        if dPmax>40
            dPmax=40;
            fprintf("==========跳出本轮==========\n")
            return
        end
        if mod(ind,100000)==0
            usetime=toc;
            fprintf("正在计算Tin=%fms，t=%.4fms，当前压强%.4f，当前最大差为%.5f，用时%f...\n",Tin,t,P,dPmax,usetime);
        end
        ind=ind+1;
    end
    fprintf("==========结束本轮==========\nTin=%fms,最大差=%.5f",Tin,dPmax)
    x=0:dt:tend;
    plot(x(1:length(p_ls)),p_ls,"g-",LineWidth=0.1)
    xlabel("时间(ms)")
    ylabel("压力(MPa)")
end
