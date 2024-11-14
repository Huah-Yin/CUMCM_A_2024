%% 初始化
% clear
clc
close all

%% 常量定义
global C A V0 get_rou get_P Aout Qout Qout2 Qin dt tend f aa0 aa1 bb1 t_of_needle h_of_needle h_dic h_max
C=0.85;
A=pi*(1.4/2)^2;
V0=pi*500*(10/2)^2;

%% 读文件
data1=xlsread("附件1-凸轮边缘曲线.xlsx");
data2=xlsread("附件2-针阀运动曲线.xlsx");
t_of_needle=[data2(:,1);data2(:,4)];
h_of_needle=[data2(:,2);data2(:,5)];

%% 拟合
f=fit(data1(:,1),data1(:,2),"fourier1");
aa0=4.826;
aa1=2.413;
bb1=-2.676e-6;

%% rou与P的关系
get_rou=@(p)3.910284764226142e-04.*p+0.808708414066112;
get_P=@(rou)(rou-0.808708414066112)/3.910284764226142e-04;

%% 流出流量函数
Aout=@(t)min((pi*(1.4/2)^2),pi*((2.5/2+hout(t)*sin(9*pi/180)*cos(9*pi/180)).^2-(2.5/2)^2)/cos(9*pi/180));
Qout=@(P1,Pout,t)(C*Aout(t)*(2*(P1-Pout)/get_rou(P1))^(0.5));
Qout2=@(P1,Pout)(C*(pi*0.7^2)*(2*(P1-Pout)/get_rou(P1))^(0.5));

%% 流入流量函数
Qin=@(P2,P1)(C*A*(2*(P2-P1)/get_rou(P2))^(0.5));

%% 转角与高度的关系
% xx=@(ttheta,ggamma)(f(ttheta).*cos(ttheta).*cos(ggamma)-f(ttheta).*sin(ttheta).*sin(ggamma));
% h_dic=[];
% for gam=0:0.01:2*pi
%     xxx=0;
%     for ttheta=1:0.01:2*pi
%         xxx_temp=xx(ttheta,gam+pi);
%         if xxx_temp>xxx
%             xxx=xxx_temp;
%         end
%     end
%     h_dic(end+1)=xxx;
%     if mod(gam,0.1)==0
%         fprintf("计算高度至%.2f/6.28\n",gam)
%     end
% end
% h_dic=h_dic-min(h_dic);
% h_max=max(h_dic);

%% 循环优化Tin
dt=0.01;   % 差分步长，单位ms
tend=10000;    % 终止时间，单位ms
t1=0.001;     % 大步长
t2=0.0001;  % 小步长
%% 大步长搜索
% p=[];
% x=0.05:t1:0.12;
% for omg=x
%     p(end+1)=get_dPmax(omg,0);
% end
% plot(x,p)
% xlabel("ω(rad/ms)")
% ylabel("Z(Mpa)")

%% 小步长搜索
% p=[];
% x=0.025:t2:0.030;
% for Tin=x
%     p(end+1)=get_dPmax(Tin,0);
% end
% plot(x,p,"b+-")
% xlabel("Tin(ms)")
% ylabel("Z(MPa)")
% [~,min_ind]=min(p);
% 
% fprintf("最优周期:%f\n",x(min_ind))

%% 灵敏性分析
V00=V0;
step=0.001*V00;
p=[];
x=V00*0.99:step:V00*1.01;
for V0=x
    p(end+1)=get_dPmax(0.0802,0);
end
plot(x,p,"r+-")
xlabel("V_{0}(mm^{3})")
ylabel("Z(MPa)")

%% 函数计算压强变化峰值
function dPmax = get_dPmax(omg,show_fig)
    global C A V0 get_rou get_P Aout Qout Qout2 Qin dt tend f aa0 aa1 bb1 t_of_needle h_of_needle h_dic h_max
    %omg-凸轮角速度
    Tin=2*pi/omg;
    Tout=50;
    P0=100; % 初始压强
    P2=0.5; % 进入压强
    m1=get_rou(P0)*V0;  % 管内实时质量，赋初值为初始质量
    m2=((f(0)-f(pi))*(pi*2.5*2.5)+20)*get_rou(0.5);
    gin=0;  % 流入转动角
    tout=0; % 流出计时器
    ind=1;
    P1=P0;   % 管内实时压强存储
    p_ls=P1;
    vin_ls=[];
    vout_ls=[];
    dPmax=0;    % 压强差最大值
    Pout=0.1;   % 外界大气压
    V1=20;
    for t=1:dt:tend
        %% 计时器重置
        if gin>=2*pi
            gin=gin-2*pi;
        end
        if tout>Tout
            tout=tout-Tout;
        end
        %% 柱塞压缩
        % 压缩高度计算
%         xx=@(ttheta,ggamma)(f(ttheta).*cos(ttheta).*cos(ggamma)-f(ttheta).*sin(ttheta).*sin(ggamma));
%         xxx=0;
%         for ttheta=1:0.01:2*pi
%             xxx_temp=xx(ttheta,gin+pi);
%             if xxx_temp>xxx
%                 xxx=xxx_temp;
%             end
%         end
%         h=xxx-f(pi);
        h=h_dic(int32(gin*100)+1);
        V2=(h_max-h)*(pi*2.5*2.5);
        if h<0.1
            %% 流入柱塞
            P2=0.5;
            m2=((f(0)-f(pi))*(pi*2.5*2.5)+20)*get_rou(0.5);
        else
            rou2=m2/(V1+V2);
            P2=get_P(rou2);
        end
        %% 流入
        if P2>P1
            vin=Qin(P2,P1);
            m2=m2-get_rou(P2)*vin*dt;
        else
            vin=0;
        end
        %% 流出
        vout=0;
        if tout<2.45
            vout=vout+Qout(P1,Pout,tout);
        end
        if vin==0 && P1>101
            vout=vout+Qout2(P1,0.5);
        end
        %% 质量更新
        if vin>0 || vout>0
            m1=m1+get_rou(P2)*vin*dt-get_rou(P1)*vout*dt;
        end
        %% 压强变化计算
        rou_now=m1/V0;
        P1=get_P(rou_now);
        if show_fig
            p_ls(end+1)=P1;
            vin_ls(end+1)=vin;
            vout_ls(end+1)=vout;
        end
        %% 计时器更新
        gin=gin+dt*omg;
        tout=tout+dt;
        %% 越界判断
        if P1<=0 || P1>=200
            dPmax=inf;
            return
        end
        if dPmax<abs(P1-P0)
            dPmax=abs(P1-P0);
        end
        if dPmax>40
            dPmax=40;
            fprintf("\n==========跳出本轮==========\n")
            return
        end
        if mod(ind,100)==0
            fprintf("\n正在计算omg=%.5frad/ms，t=%.4fms，当前压强P1=%.4fMPa;P2=%.4fMpa;压缩率%.4f;Vin=%.4f;Vout=%.4f;h=%.4f,V2=%.4f，当前最大差为%.5f...",omg,t,P1,P2,gin/pi,vin,vout,h,V2,dPmax);
        end
        ind=ind+1;
    end
    fprintf("\n==========结束本轮==========\nomg=%fms,最大差=%.5f",Tin,dPmax)
    if show_fig
        fprintf("\n方差为%f\n",var(p_ls))
        range=40000;
        figure
        x=0:dt:tend;
        plot(x(1:range),p_ls(1:range),"g-",LineWidth=0.1)
        xlabel("时间(ms)")
        ylabel("压力(MPa)")

        figure
        plot(x(1:length(p_ls)),p_ls,"g-",LineWidth=0.1)
        xlabel("时间(ms)")
        ylabel("压力(MPa)")

        figure
        plot(x(1:range),vin_ls(1:range),"g-",x(1:range),vout_ls(1:range),"r-");
        xlabel("时间(ms)")
        ylabel("流量(mm^{3}/ms)")
        legend("流入流量","流出流量")
    end
end

%% 针阀升程与时间的关系
function h=hout(t)
    global t_of_needle h_of_needle
    if t>=0.45-1e-4 && t<=2+1e-4
        h=2;
        return
    elseif t>=2.46-1e-4
        h=0;
        return
    else
        h=h_of_needle(find(abs(t_of_needle-t)<1e-6));
        return
    end
end

