%% 初始化
clear
clc
close all
%% 参数赋值
global data chang kuan h2 h1 as gs d h1
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
yita_cos_day=zeros(1,12);
yita_at_day=zeros(1,12);
yita_trunc_day=zeros(1,12);
yita_ref_day=ones(1,12).*0.92;
%% 各镜数据存储
len=length(data);
res_of_all=zeros(len,3);
res_of_all(:,1:2)=data;
rec_time=0;
%% 反射光锥参量计算
theta_dot=atan(R1/d)*180/pi;
%% 遍历时间
D_ind=0;
for D=[-59,-28,0,31,61,92,122,153,184,214,245,275]
    D_ind=D_ind+1;
    ST_ind=0;
    %% 结果初始化
    yita_sb_time=zeros(1,5);
    yita_cos_time=zeros(1,5);
    yita_at_time=zeros(1,5);
    yita_trunc_time=zeros(1,5);
    E_time=zeros(1,5);
%     for ST=[9]
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
        yita_cos=zeros(len,1);
        yita_at=zeros(len,1);
        yita_trunc=ones(len,1).*(-1);
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
                    %% 截断效率计算
                    if black==1
                        for line_ind=1:5
                            det=4*(x1*s{line_ind}(1)+y1*s{line_ind}(2))^2-4*(s{line_ind}(1)^2+s{line_ind}(2)^2)*(x1^2+y1^2-R0^2);
                            if det<0
                                unin_num=unin_num+1;    % 不能射入
                            else
                                t_sol=(-2*(x1*s{line_ind}(1)+y1*s{line_ind}(2))-sqrt(det))/(2*(s{line_ind}(1)^2+s{line_ind}(2)^2));
                                z_sol=z1+s{line_ind}(3)*t_sol;
                                
                                if z_sol>=h2-h3/2 && z_sol<=h2+h3/2
                                    in_num=in_num+1;    % 能射入
                                    in{num}(ind)=0;
                                else
                                    unin_num=unin_num+1;    % 不能射入
                                    in{num}(ind)=1;
                                end
                            end
                        end
                    end
                    %% 计算余弦效率
                    n_f=[cos(am*pi/180)*sin(gm*pi/180);cos(am*pi/180)*cos(gm*pi/180);sin(am*pi/180)];
                    cosn(ind)=real(dot(s_f,n_f));
                    ind=ind+1;
                end
            end
            %% 统计阴影遮挡效率
            yita_sb(num)=1-bla_num/(bla_num+unbla_num);
            %% 计算单板平均余弦效率
            yita_cos(num)=mean(cosn);
            %% 计算大气透射率
            dHR=sqrt(x^2+y^2+(h2-h1)^2);
            yita_at(num)=0.99321-0.0001176*dHR+1.97e-8*dHR^2;
            %% 计算截断效率
            if yita_sb(num)~=0
                yita_trunc(num)=1-in_num/(in_num+unin_num);
            else
                yita_trunc(num)=0;
            end
            toc
            yita=yita_sb(num)*yita_cos(num)*yita_at(num)*yita_trunc(num)*yita_ref_day(D_ind);
            res_of_all(num,3)=yita;
            E=E+chang*kuan*yita;
            fprintf('完成计算%d月%.1f时第%d块定日镜...\n',D_ind,ST,num);
        end
        E=E*DNI;
        yita_sb_time(ST_ind)=mean(yita_sb);
        yita_cos_time(ST_ind)=mean(yita_cos);
        yita_at_time(ST_ind)=mean(yita_at);
        yita_trunc_time(ST_ind)=mean(yita_trunc);
        E_time(ST_ind)=E/(chang*kuan*len);
        rec_time=rec_time+1;
    end
    yita_sb_day(D_ind)=mean(yita_sb_time);
    yita_cos_day(D_ind)=mean(yita_cos_time);
    yita_at_day(D_ind)=mean(yita_at_time);
    yita_trunc_day(D_ind)=mean(yita_trunc_time);
    E_day=mean(E_time);
    %% 结果整理
    res{D_ind,1}=yita_sb_day(D_ind)*yita_cos_day(D_ind)*yita_at_day(D_ind)*yita_trunc_day(D_ind)*yita_ref_day(D_ind);
    res{D_ind,2}=yita_cos_day(D_ind);
    res{D_ind,3}=yita_sb_day(D_ind);
    res{D_ind,4}=yita_trunc_day(D_ind);
    res{D_ind,5}=E_day;
    fprintf("数据已记录至res与res_of_all\n")
    xlswrite("2023A-定日镜问题1结论.xlsx",res_of_all)
    fprintf("文件已保存\n")
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
    global data chang kuan h1
    %% 寻找最近定日镜
    distance=((data(:,1)-data(now_ind,1)).^2 + (data(:,2)-data(now_ind,2)).^2).^0.5;
    [~,near_ind]=mink(distance,10);
%     near_ind(find(distance(near_ind)<0.001))=[];  % 除去自身
%     near_ind=1:length(data);
%     near_ind=near_ind';
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
%             fprintf("xb=%f,yb=%f\n",xb,yb)
            if ((-kuan/2<=xb && xb<=kuan/2) && (-chang/2<=yb && yb<=chang/2))
                res=0;
                fprintf("挡")
                return
            end
%%
%             D=[cos(am)*sin(gm)+s_now(2)*cos(gm)/s_now(1), sin(am);
%                sin(am)*sin(gm)+s_now(3)*cos(gm)/s_now(1), -cos(am)];
%             D1=[(M(1)+x0)*s_now(2)/s_now(1)-M(2)+y0, sin(am);
%                 (M(1)+x0)*s_now(3)/s_now(1)-M(3)+h1, -cos(am)];
%             D2=[cos(am)*sin(gm)+s_now(2)*cos(gm)/s_now(1), (M(1)+x)*s_now(2)/s_now(1)-M(2);
%                sin(am)*sin(gm)+s_now(3)*cos(gm)/s_now(1),  (M(1)+h1)*s_now(3)/s_now(1)-M(3)];
%             dD1=det(D1);
%             dD2=det(D2);
%             dD=det(D);
%             a=dD1/dD;
%             b=dD2/dD;
            %% 测试
%             A=[cos(gm*pi/180),sin(gm*pi/180),0;-sin(gm*pi/180),cos(gm*pi/180),0;0,0,1];
%             B=[1,0,0;0,cos(am*pi/180),-sin(am*pi/180);0,sin(am*pi/180),cos(am*pi/180)];
%             p1=B*A*[3;0;3];
%             p2=B*A*[3;0;-3];
%             p3=B*A*[-3;0;-3];
%             p4=B*A*[-3;0;3];
%             t=a*(cos(gm*pi/180)-M(1))/s_now(1);
%             fill3([p1(1)+x,p2(1)+x,p3(1)+x,p4(1)+x,p1(1)+x],[p1(2)+h1,p2(2)+h1,p3(2)+h1,p4(2)+h1,p1(2)+h1],[p1(3)+y,p2(3)+y,p3(3)+y,p4(3)+y,p1(3)+y],"r",[M(1)+(t-10)*s_now(1),M(1)+t*s_now(1)],[M(2)+(t-10)*s_now(2),M(2)+t*s_now(2)],[M(3)+(t-10)*s_now(3),M(3)+t*s_now(3)],"b")
%             grid on
%             fprintf("a=%f;b=%f\n",a,b);
            %% 结束
            %%
            if (a<chang/2 && a>-chang/2) && (b<kuan/2 && b>-kuan/2)
                fprintf("a=%f;b=%f\n",a,b);
                fprintf("=================>遮挡")
                res=0;
                return
            end
        end
    end
    res=1;
    return
end
