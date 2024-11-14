%% 初始化
clear
clc
close all

%% 主函数调用
global result sen_idx sb_num_sen
x=9:13;
result=zeros(2,length(x));
sen_idx=1;
for i=x
    sb_num_sen=i;
    main();
    sen_idx=sen_idx+1;
end
fprintf("结果已存至变量result中")

%% ======主函数======
function main()
    % sb_num:考虑遮挡的定日镜数量
    %% 参数赋值
    global data chang kuan h2 h1 as gs d h1 sb_num_sen
    global result sen_idx
    step=0.5;   % 网格生成步长
    H=3;
    h1=4;
    h2=80;
    h3=8;
    R0=3.5;
    Rm=100;
    chang=6;
    kuan=6;
    d=149597870700; % 日地距离
    R1=6.963e8; % 太阳半径
    fai=39.4;
    q=0.05;   % 遮挡角度允许误差（角度）
    %% 导入数据
    data=xlsread('附件.xlsx');
    res={'平均光学效率','平均余弦效率','平均阴影遮挡效率','平均截断效率','单位平均热功率'};
    %% 结果初始化
    yita_sb_day=zeros(1,12);
    %% 各镜数据存储
    len=length(data);
    res_of_all=zeros(len,3);
    res_of_all(:,1:2)=data;
    rec_time=0;
    %% 反射光锥参量计算
    theta_dot=atan(R1/d)*180/pi;
    %% 遍历时间
    D_ind=0;
    for D=[0,184]
        D_ind=D_ind+1;
        ST_ind=0;
        %% 结果初始化
        yita_sb_time=zeros(1,5);
        for ST=[9,10.5,12,13.5,15]
            ST_ind=ST_ind+1;
            E=0;    % 输出热功率初始化
            %% 计算太阳高度角与方位角
            de=asin(sin(2*pi*D/350)*sin(2*pi*23.45/360));
            omg=pi*(ST-12)/12;
            as=asin(cos(de)*cos(fai*pi/180)*cos(omg)+sin(de)*sin(fai*2*pi/360))*180/pi;          % 太阳高度角
            gs=acos((sin(de)-sin(as*pi/180)*sin(fai*2*pi/360))/cos(as*pi/180)*cos(fai*2*pi/360))*180/pi;  % 太阳方位角
            %% 计算DNI
            DNI=1.366*((0.4237-0.00821*(6-H)^2)+(0.5055+0.00595*(6.5-H)^2)*exp(-(0.2711+0.01858*(2.5-H)^2)/sin(as*pi/180)));
            %% 遍历定日镜
            in=cell(1,len);
            yita_sb=zeros(len,1);
            for num=1:len
                tic
                x=data(num,1);
                y=data(num,2);
                [am,gm]=get_amgm(x,y);
                bla_num=0;
                unbla_num=0;
                in_num=0;
                unin_num=0;
                %% 反射光锥轴线向量计算
                s0=-[cos(as*pi/180)*sin(gs*pi/180);cos(as*pi/180)*cos(gs*pi/180);sin(as*pi/180)];
                n=[cos(am*pi/180)*sin(gm*pi/180);cos(am*pi/180)*cos(gm*pi/180);sin(am*pi/180)];
                s_dot=s0+2*abs(dot(s0,n))*n;
                s_dot=s_dot./norm(s_dot);
                s{5}=s_dot;
                omg_dot=asin(s_dot(3))*180/pi;
                %% 上下采样线方向向量计算
                s{1}=s_dot+[0;0;sin(theta_dot*pi/180)/sin((90-omg_dot-theta_dot)*pi/180)];
                s{1}=s{1}./norm(s{1});   % 单位化
                s{2}=s_dot+[0;0;sin(theta_dot*pi/180)/sin((90+omg_dot-theta_dot)*pi/180)];
                s{2}=s{2}./norm(s{2});
                %% 左右采样线方向向量计算
                m1=s{1}./cos((theta_dot+omg_dot)*pi/180)-s_dot;
                temp_cr=cross(s_dot,m1);
                m4=temp_cr.*norm(m1)./norm(temp_cr);
                m3=-m4;
                s{3}=m3+s_dot;
                s{3}=s{3}./norm(s{3});
                s{4}=m4+s_dot;
                s{4}=s{4}./norm(s{4});
                %% 划分网格
                ind_max=floor(chang/step+1)*floor(kuan/step+1);
                cosn=zeros(1,ind_max);
                ind=1;
                for i=-chang/2:step:chang/2
                    for k=-kuan/2:step:kuan/2
                        x0=i;
                        y0=0;
                        z0=k;
                        %% 变换
                        A=[cos(gm*pi/180),sin(gm*pi/180),0;-sin(gm*pi/180),cos(gm*pi/180),0;0,0,1];
                        B=[1,0,0;0,cos(am*pi/180),-sin(am*pi/180);0,sin(am*pi/180),cos(am*pi/180)];
                        t=B*A*[x0;0;z0];
                        x1=t(1)+x;
                        y1=t(2)+y;
                        z1=t(3)+h1;
                        M=[x1;y1;z1];
                        %% 吸收塔遮挡判断
                        s_f=[cos(as*pi/180)*sin(gs*pi/180);cos(as*pi/180)*cos(gs*pi/180);sin(as*pi/180)];
                        if x1<0 && y1>0   % 仅判断第二象限是否被遮挡
                            jijiao=(atan(y1/x1)+pi)*180/pi; % 极角
                            theta_s=(180-gs)+90;
                            det_j=atan(R0/Rm)*180/pi+q;
                            if jijiao>theta_s-det_j && jijiao<theta_s+det_j
                                det=4*(x1*s_f(1)+y1*s_f(2))^2-4*(s_f(1)^2+s_f(2)^2)*(x1^2+y1^2-R0^2);
                                if det<0
                                    %% 定日镜遮挡判断
                                    if jud_sb(s_f,s,num,M)
                                        black=1;  % 不被遮挡
                                        unbla_num=unbla_num+1;
                                    else
                                        black=0;  % 被遮挡
                                        bla_num=bla_num+1;
                                    end
                                else
                                    t_sol=(-2*(x1*s_f(1)+y1*s_f(2))-sqrt(det))/(2*(s_f(1)^2+s_f(2)^2));
                                    z_sol=z1+s_f(3)*t_sol;
                                    if z_sol>=0 && z_sol<=h2+h3/2
                                        black=0;  % 被遮挡
                                        bla_num=bla_num+1;
                                    else
                                        %% 定日镜遮挡判断
                                        if jud_sb(s_f,s,num,M)
                                            black=1;  % 不被遮挡
                                            unbla_num=unbla_num+1;
                                        else
                                            black=0;  % 被遮挡
                                            bla_num=bla_num+1;
                                        end
                                    end
                                end
                            else
                                %% 定日镜遮挡判断
                                if jud_sb(s_f,s,num,M)
                                    black=1;  % 不被遮挡
                                    unbla_num=unbla_num+1;
                                else
                                    black=0;  % 被遮挡
                                    bla_num=bla_num+1;
                                end
                            end
                        else
                            %% 定日镜遮挡判断
                            if jud_sb(s_f,s,num,M)
                                black=1;  % 不被遮挡
                                unbla_num=unbla_num+1;
                            else
                                black=0;  % 被遮挡
                                bla_num=bla_num+1;
                            end
                        end
                        ind=ind+1;
                    end
                end
                %% 统计阴影遮挡效率
                yita_sb(num)=1-bla_num/(bla_num+unbla_num);
                toc
                fprintf('完成计算%d月%.1f时第%d块定日镜，遮挡参数%d...\n',D_ind,ST,num,sb_num_sen);
            end
            yita_sb_time(ST_ind)=mean(yita_sb);
            rec_time=rec_time+1;
        end
        result(D_ind,sen_idx)=mean(yita_sb_time);
    end
end

%% 计算定日镜高度角、方向角
function [am,gm]=get_amgm(x,y)
    global h2 h1 as gs d
    r=sqrt(x^2+y^2);
    %% 计算定日镜高度角
    theta1=atan((h2-h1)/r)*180/pi;
    am=theta1+(as-theta1)/2;
    %% 计算定日镜方位角
    if y>0 && ( x>=0 || -atan(x/y)*180/pi+gs<=180 ) % 情况1
        beta=180*atan(y/x)/pi;
        theta2=beta+gs-sign(beta)*90;
        l=sqrt(r^2+d^2-2*r*d*cos(theta2*pi/180));
        theta3=180*acos((r^2+l^2-d^2)/(2*r*l))/pi;
        if x>=0
            gm=180-(theta3/2-(90-beta));
        else
            gm=-(theta3/2-(90-beta));
        end
    elseif y<0 && ( x<=0 || -atan(x/y)*180/pi+gs<=180 ) % 情况2
        beta=180*atan(x/y)/pi;
        theta4=beta+180-gs;
        l=sqrt(r^2+d^2-2*r*d*cos(theta4*pi/180));
        gm=beta+1/2*180*acos((r^2+l^2-d^2)/(2*r*l))/pi;
    elseif  y>=0 && ( x<=0 && -atan(x/y)*180/pi+gs>=180 ) % 情况3
        beta=180*abs(atan(x/y))/pi;
        theta5=(180-beta)+(180-gs);
        l=sqrt(r^2+d^2-2*r*d*cos(theta5*pi/180));
        gm=180-(beta-1/2*180/pi*acos((r^2+l^2-d^2)/(2*r*l)));
    else % 情况4
        beta=180*atan(x/y)/pi;
        theta6=180-gs+beta;
        l=sqrt(r^2+d^2-2*r*d*cos(theta6*pi/180));
        gm=180-1/2*180/pi*acos((r^2+l^2-d^2)/(2*r*l))+beta;
    end
end

%% 判断是否被定日镜遮挡：1-不被遮挡，0-被遮挡
function res=jud_sb(s_f,s,now_ind,M)
    global data chang kuan h1 sb_num_sen
    %% 寻找最近定日镜
    distance=((data(:,1)-data(now_ind,1)).^2 + (data(:,2)-data(now_ind,2)).^2).^0.5;
    [~,near_ind]=mink(distance,sb_num_sen);
    %% 遍历遮挡定日镜
    for line_ind=0:5
        if line_ind==0
            s_now=s_f;
        else
            s_now=s{line_ind};
        end
        for ind=near_ind'
            if ind==now_ind
                continue
            end
            x0=data(now_ind,1);
            y0=data(now_ind,2);
            [am0,gm0]=get_amgm(x0,y0);
            x=data(ind,1);
            y=data(ind,2);
            [am,gm]=get_amgm(x,y);
            A0=[cos(gm0*pi/180),sin(gm0*pi/180),0;-sin(gm0*pi/180),cos(gm0*pi/180),0;0,0,1];
            B0=[1,0,0;0,cos(am0*pi/180),-sin(am0*pi/180);0,sin(am0*pi/180),cos(am0*pi/180)];
            T0=B0*A0;
%             T0i=inv(T0);
            A=[cos(gm*pi/180),sin(gm*pi/180),0;-sin(gm*pi/180),cos(gm*pi/180),0;0,0,1];
            B=[1,0,0;0,cos(am*pi/180),-sin(am*pi/180);0,sin(am*pi/180),cos(am*pi/180)];
            T=B*A;
%             og=T0*M+[x0;y0;h1];
            b=T0\(M-[x;y;h1]);
            ss=T0\s_now;
            xb=(ss(3)*b(1)-ss(1)*b(3))/ss(3);
            yb=(ss(3)*b(2)-ss(2)*b(3))/ss(3);
            if ((-kuan/2<=xb && xb<=kuan/2) && (-chang/2<=yb && yb<=chang/2))
                res=0;
                fprintf("挡")
                return
            end
        end
    end
    res=1;
    return
end
