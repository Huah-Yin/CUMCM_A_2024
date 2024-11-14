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
% T_best=[];
% Z_in_best=[];
% for tc=[5,10]
%     p=[];
%     x=10.2:t1:11;
%     for Tin1=x
%         p(end+1)=get_dPmax(Tin1,10.7602,tc*1000,0);
%     end
%     figure
%     plot(x(1:length(p)),p)
%     title(strcat(num2str(tc),"s"))
%     xlabel("Tin1(ms)")
%     ylabel("Z(Mpa)")
%     [Z_in_best(end+1),T_best_ind]=min(p);
%     T_best(end+1)=x(T_best_ind);
% end

%% 小步长搜索
% T_best=[];
% Z_in_best=[];
% for tc=[5,10]
%     p=[];
%     x=10.72:t2:10.78;
%     for Tin1=x
%         p(end+1)=get_dPmax(Tin1,10.7602,tc*1000,0);
%     end
%     figure
%     plot(x(1:length(p)),p)
%     title(strcat(num2str(tc),"s"))
%     xlabel("Tin1(ms)")
%     ylabel("Z(Mpa)")
%     [Z_in_best(end+1),T_best_ind]=min(p);
%     T_best(end+1)=x(T_best_ind);
% end

get_dPmax(10.9568,10.7602,2*1000,1)
get_dPmax(10.7570,10.7602,5*1000,1)
get_dPmax(10.7434,10.7602,10*1000,1)


%% 报告打印
% fprintf("==========计算结束==========\n")
% tls=[5,10];
% for i=1:3
%     fprintf("%d秒:Tin1=%.6f,Z=%.4f\n",tls(i),T_best(i),Z_in_best(i))
% end
% figure
% plot(tls,T_best,"k-o")
% title("Tin与调压时间的关系图")
% xlabel("调压时间")
% ylabel("Tin1(ms)")
% xticks(tls)
% fprintf("===========================\n")
%% 函数计算压强变化峰值
function dPmax = get_dPmax(Tin1,Tin2,t_change,show_fig)
    global C A V0 get_rou get_P Qout Qin dt tend
    %Tin-流入周期
    Tout=100;
    P0=100; % 初始压强
    Pw=150; % 目标压强
    P2=160; % 进入压强
    m=get_rou(P0)*V0;  % 实时质量，赋初值为初始质量
    tin=1;  % 流入计时器
    tout=1; % 流出计时器
    ind=1;
    P=P0;   % 实时压强存储
    if show_fig
        p_ls=P0;
    end
    dPmax=0;    % 压强差最大值
    global data_out
    data_out={};
    for t=1:dt:tend
        tic
        %% 周期确定
        if t>t_change
            Tin=Tin2;
        else
            Tin=Tin1;
        end
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
        if show_fig && mod(ind,1000)==0
            p_ls(end+1)=P;
        end
        %% 计时器更新
        tin=tin+dt;
        tout=tout+dt;
        %% 越界判断
        if P<=0 || P>=200
            dPmax=inf;
            fprintf("==========跳出本轮（越界）==========\n")
            return
        end
        if t>t_change
            if dPmax<abs(P-Pw)
                dPmax=abs(P-Pw);
            end
        end
        if dPmax>40
            dPmax=40;
            fprintf("==========跳出本轮==========\n")
            return
        end
        if mod(ind,100000)==0
            usetime=toc;
            fprintf("正在计算%d秒：Tin1=%fms，t=%.4fms，当前压强%.4f，当前最大差为%.5f，用时%f...\n",t_change/1000,Tin1,t,P,dPmax,usetime);
        end
        ind=ind+1;
    end
    fprintf("==========结束本轮==========\nTin1=%fms,最大差=%.5f",Tin1,dPmax)
    if show_fig
        x=0:dt:tend;
        figure
        plot(x(1:length(p_ls)),p_ls,"g-",LineWidth=0.1)
        xlabel("时间(ms)")
        ylabel("压力(MPa)")
        title(strcat(num2str(t_change/1000),"s"))
    end
end
